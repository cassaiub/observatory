import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import build_dq, write_mef
from cassa_photometry.phase4_diagnostics import pipeline, psf, stages

_FWHM = 4.0
_SIGMA = _FWHM / 2.3548200450309493


def _star_image(shape=(200, 200), n=18, fwhm=_FWHM, flux=6000.0, noise=5.0, seed=1):
    """Render isolated Gaussian stars of known FWHM on a noisy background."""
    rng = np.random.default_rng(seed)
    img = rng.normal(100.0, noise, size=shape).astype(np.float32)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    sigma = fwhm / 2.3548200450309493
    positions = []
    margin = 20
    for _ in range(n):
        x = rng.uniform(margin, shape[1] - margin)
        y = rng.uniform(margin, shape[0] - margin)
        if any(np.hypot(x - px, y - py) < 30 for px, py in positions):
            continue
        img += flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sigma ** 2))
        positions.append((x, y))
    return img, positions


def _make_wcs_header(shape):
    w = WCS(naxis=2)
    w.wcs.crpix = [shape[1] / 2, shape[0] / 2]
    w.wcs.crval = [339.27, 34.42]
    w.wcs.cdelt = [-2.0 / 3600, 2.0 / 3600]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    h = w.to_header()
    h["PIXSCALE"] = 2.0
    h["STACKCNT"] = 4
    h["FILTER"] = "R"
    return h


def test_estimate_fwhm_recovers_known():
    img, _ = _star_image()
    res = psf.estimate_fwhm(img, fwhm_guess=3.5, threshold=5.0)
    assert res["n_used"] >= 5
    assert abs(res["fwhm_px"] - _FWHM) < 0.7   # recovered within ~0.7 px
    assert res["radii"] is not None and res["profile"] is not None


def test_diagnose_raw(tmp_path):
    img, _ = _star_image()
    path = tmp_path / "raw.fits"
    h = fits.Header()
    h["OBJECT"] = "SYN"
    h["FILTER"] = "R"
    fits.PrimaryHDU(img, h).writeto(str(path))
    metrics, figs = stages.diagnose_raw(str(path), load_config(), _logger())
    assert np.isfinite(metrics["fwhm_px"]) and len(figs) == 1


def test_diagnose_calibrated(tmp_path):
    img, _ = _star_image()
    err = np.sqrt(np.clip(img, 0, None)).astype(np.float32)
    dq = build_dq(img.shape, saturated=(img > img.max() * 0.99))
    path = tmp_path / "calibrated_syn.fits"
    h = fits.Header()
    h["FILTER"] = "R"
    h["BUNIT"] = "electron"
    write_mef(str(path), sci=img, err=err, dq=dq, header=h)
    metrics, figs = stages.diagnose_calibrated(str(path), load_config(), _logger())
    assert "median_err" in metrics and len(figs) == 1


def test_diagnose_master_with_wcs(tmp_path):
    img, _ = _star_image()
    err = np.sqrt(np.clip(img, 0, None)).astype(np.float32)
    path = tmp_path / "Master_syn.fits"
    write_mef(str(path), sci=img, err=err, dq=build_dq(img.shape),
              header=_make_wcs_header(img.shape))
    metrics, figs = stages.diagnose_master(str(path), load_config(), _logger())
    assert metrics["wcs_solved"] is True and len(figs) == 1


def test_diagnose_photometry(tmp_path):
    rng = np.random.default_rng(0)
    n = 200
    mag = rng.uniform(12, 21, n)
    magerr = 0.01 * 10 ** (0.4 * (mag - 15))
    tbl = Table({
        "NUMBER": np.arange(n), "ALPHA_J2000": rng.uniform(339, 339.1, n),
        "DELTA_J2000": rng.uniform(34.3, 34.5, n),
        "X_IMAGE": rng.uniform(0, 1024, n), "Y_IMAGE": rng.uniform(0, 1024, n),
        "FLUX_ISO": 10 ** (-0.4 * (mag - 25)),
        "FLUXERR_ISO": rng.uniform(1, 50, n),
        "MAG_ISO": mag, "MAGERR_ISO": magerr,
        "SNR": 1.0857 / magerr, "ISOAREA_IMAGE": rng.uniform(5, 200, n),
        "ELLIPTICITY": rng.uniform(0, 0.4, n), "FLAGS": np.zeros(n, int),
    })
    cat = tmp_path / "syn_catalog.csv"
    tbl.write(str(cat), format="csv", overwrite=True)
    metrics, figs = stages.diagnose_photometry(str(cat), load_config(), _logger())
    assert metrics["n_sources"] == n and len(figs) == 1


def test_full_report(tmp_path):
    """Assemble a mini work dir (phase1/phase2) and confirm the report is produced."""
    img, _ = _star_image()
    err = np.sqrt(np.clip(img, 0, None)).astype(np.float32)
    p1_dir = tmp_path / "phase1"
    p2_dir = tmp_path / "phase2"
    p1_dir.mkdir()
    p2_dir.mkdir()

    fits.PrimaryHDU(img, fits.Header({"FILTER": "R", "OBJECT": "SYN"})).writeto(str(tmp_path / "raw_R.fits"))
    write_mef(str(p1_dir / "calibrated_R.fits"), sci=img, err=err,
              dq=build_dq(img.shape), header=fits.Header({"FILTER": "R", "BUNIT": "electron"}))
    write_mef(str(p2_dir / "Master_R.fits"), sci=img, err=err,
              dq=build_dq(img.shape), header=_make_wcs_header(img.shape))

    metrics = pipeline.run(run_dir=str(p2_dir), raw=str(tmp_path / "raw_R.fits"),
                           config=load_config())
    diag = tmp_path / "phase4"  # sibling of phase2, written directly (no nested dir)
    assert (diag / "diagnostics_report.pdf").exists()
    assert (diag / "metrics.json").exists()
    assert "stage_2_master" in metrics
    # the calibrated frame is found in the sibling phase1 directory
    assert metrics.get("stage_1_calibrated") is not None


def test_full_report_legacy_layout(tmp_path):
    """The older nested layout (calibrated frames in the run dir's parent) still resolves."""
    img, _ = _star_image()
    err = np.sqrt(np.clip(img, 0, None)).astype(np.float32)
    run_dir = tmp_path / "Run01"
    run_dir.mkdir()

    fits.PrimaryHDU(img, fits.Header({"FILTER": "R", "OBJECT": "SYN"})).writeto(str(tmp_path / "raw_R.fits"))
    write_mef(str(tmp_path / "calibrated_R.fits"), sci=img, err=err,
              dq=build_dq(img.shape), header=fits.Header({"FILTER": "R", "BUNIT": "electron"}))
    write_mef(str(run_dir / "Master_R.fits"), sci=img, err=err,
              dq=build_dq(img.shape), header=_make_wcs_header(img.shape))

    metrics = pipeline.run(run_dir=str(run_dir), raw=str(tmp_path / "raw_R.fits"),
                           outdir=str(run_dir / "diagnostics"), config=load_config())
    assert (run_dir / "diagnostics" / "diagnostics_report.pdf").exists()
    assert metrics.get("stage_1_calibrated") is not None


def _logger():
    from cassa_photometry.logging_utils import get_logger
    return get_logger("test_phase4")
