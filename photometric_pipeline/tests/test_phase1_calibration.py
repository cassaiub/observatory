import glob
import os

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from conftest import FixtureProfile

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import DQ_COSMIC_RAY
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.phase1_calibration.data_models import load_standardized_ccds
from cassa_photometry.phase1_calibration.pipeline import (
    _build_bpm,
    _outlier_threshold,
    _stack_scatter,
    build_master_bias,
    build_master_dark,
    build_master_flat,
)
from cassa_photometry.phase1_calibration.pipeline import (
    run as run_phase1,
)
from cassa_photometry.phase1_calibration.processor import (
    UniversalProcessor,
    crop_master_to_frame,
    subtract_scaled_dark,
)

LOGGER = get_logger("test_phase1")
SHAPE = (16, 16)


def _write(tmp_path, name, value, img_type, exptime, filt="R"):
    header = fits.Header()
    header["IMAGETYP"] = img_type
    header["EXPTIME"] = exptime
    header["EGAIN"] = 1.0
    header["READNOIS"] = 5.0
    header["FILTER"] = filt
    header["TELESCOP"] = "T32"
    data = np.full(SHAPE, value, dtype=np.float32)
    path = tmp_path / name
    fits.PrimaryHDU(data=data, header=header).writeto(str(path))
    return str(path)


def _master(data, exptime=None, gain=None, bias_subtracted=None, ncombine=None, scatter=None):
    """A stand-in master frame carrying the provenance the pipeline records."""
    ccd = CCDData(np.asarray(data, dtype=np.float64), unit=u.adu)
    if exptime is not None:
        ccd.meta["EXPTIME"] = float(exptime)
    if gain is not None:
        ccd.meta["GAINVAL"] = float(gain)
    if bias_subtracted is not None:
        ccd.meta["BIASSUB"] = bool(bias_subtracted)
    if ncombine is not None:
        ccd.meta["NCOMBINE"] = int(ncombine)
    if scatter is not None:
        ccd.uncertainty = StdDevUncertainty(np.asarray(scatter) / np.sqrt(ncombine))
    return ccd


# --- The dark-exposure bug ----------------------------------------------------

def test_dark_is_rescaled_to_the_science_exposure():
    """A 120s master dark applied to a 60s frame must be halved, not subtracted whole."""
    science = CCDData(np.full(SHAPE, 1000.0), unit=u.adu)
    master_dark = _master(np.full(SHAPE, 100.0), exptime=120.0, bias_subtracted=True)

    result = subtract_scaled_dark(science, master_dark, data_exposure=60.0, logger=LOGGER)

    assert np.allclose(result.data, 950.0)      # 1000 - 100*(60/120)
    assert not np.allclose(result.data, 900.0)  # the unscaled result this used to give


def test_dark_is_not_rescaled_when_exposures_match():
    science = CCDData(np.full(SHAPE, 1000.0), unit=u.adu)
    master_dark = _master(np.full(SHAPE, 100.0), exptime=60.0, bias_subtracted=True)

    result = subtract_scaled_dark(science, master_dark, data_exposure=60.0, logger=LOGGER)

    assert np.allclose(result.data, 900.0)


def test_scaling_an_unsubtracted_dark_warns(caplog):
    """Scaling a dark that still holds the bias pedestal must not happen silently."""
    science = CCDData(np.full(SHAPE, 1000.0), unit=u.adu)
    master_dark = _master(np.full(SHAPE, 100.0), exptime=120.0, bias_subtracted=False)

    with caplog.at_level("WARNING"):
        subtract_scaled_dark(science, master_dark, data_exposure=60.0, logger=LOGGER)

    assert any("not bias-subtracted" in r.message for r in caplog.records)


def test_master_dark_without_exposure_falls_back_to_unscaled(caplog):
    science = CCDData(np.full(SHAPE, 1000.0), unit=u.adu)
    master_dark = _master(np.full(SHAPE, 100.0))  # no EXPTIME recorded

    with caplog.at_level("WARNING"):
        result = subtract_scaled_dark(science, master_dark, data_exposure=60.0, logger=LOGGER)

    assert np.allclose(result.data, 900.0)
    assert any("no usable EXPTIME" in r.message for r in caplog.records)


def test_processor_rescales_a_mismatched_master_dark(tmp_path):
    """End-to-end guard: the exposure mismatch used to be invisible inside ISR."""
    from cassa_photometry.phase1_calibration.processor import UniversalProcessor

    cfg = load_config()
    instrument = FixtureProfile()
    path = _write(tmp_path, "sci.fits", 1000.0, "Light", 60.0)
    std_ccds = load_standardized_ccds(path, instrument, cfg)

    # 120s master dark holding 100 ADU; the 60s science frame should lose 50.
    master_dark = _master(np.full(SHAPE, 100.0), exptime=120.0, bias_subtracted=True)
    processor = UniversalProcessor(master_dark=master_dark, config=cfg, logger=LOGGER)

    frame = processor.process_science_frame(std_ccds)[0]

    assert np.allclose(np.median(frame.ccd.data), 950.0)


# --- Master-frame provenance --------------------------------------------------

def test_master_dark_records_exposure_and_bias_provenance(tmp_path):
    cfg = load_config()
    instrument = FixtureProfile()
    biases = [_write(tmp_path, f"b{i}.fits", 100.0, "Bias", 0.0) for i in range(3)]
    darks = [_write(tmp_path, f"d{i}.fits", 130.0, "Dark", 60.0) for i in range(3)]

    m_bias = build_master_bias(biases, instrument, cfg, LOGGER)
    m_dark = build_master_dark(darks, instrument, m_bias, cfg, LOGGER)

    assert m_dark.meta["EXPTIME"] == 60.0
    assert m_dark.meta["BIASSUB"] is True
    assert m_dark.meta["NCOMBINE"] == 3
    assert np.allclose(m_dark.data, 30.0)  # bias pedestal removed


def test_mixed_dark_exposures_are_normalised_before_combining(tmp_path):
    """Darks of 60s and 120s must combine into one self-consistent master."""
    cfg = load_config()
    instrument = FixtureProfile()
    biases = [_write(tmp_path, f"b{i}.fits", 100.0, "Bias", 0.0) for i in range(3)]
    m_bias = build_master_bias(biases, instrument, cfg, LOGGER)

    # 0.5 ADU/s of dark current, sampled at two exposure times.
    darks = [_write(tmp_path, f"d60_{i}.fits", 130.0, "Dark", 60.0) for i in range(2)]
    darks += [_write(tmp_path, f"d120_{i}.fits", 160.0, "Dark", 120.0) for i in range(2)]

    m_dark = build_master_dark(darks, instrument, m_bias, cfg, LOGGER)

    ref = m_dark.meta["EXPTIME"]
    assert np.allclose(m_dark.data, 0.5 * ref), "rescaled frames should agree on one rate"


def test_master_flat_is_dark_subtracted(tmp_path):
    """Dark current in the flats must not be normalised into the flat field."""
    cfg = load_config()
    instrument = FixtureProfile()
    biases = [_write(tmp_path, f"b{i}.fits", 100.0, "Bias", 0.0) for i in range(3)]
    m_bias = build_master_bias(biases, instrument, cfg, LOGGER)
    m_dark = _master(np.full(SHAPE, 600.0), exptime=60.0, gain=1.0, bias_subtracted=True)

    flats = [_write(tmp_path, f"f{i}.fits", 1100.0, "Flat", 10.0) for i in range(3)]
    m_flat = build_master_flat(flats, instrument, m_bias, m_dark, cfg, LOGGER)

    # 1100 - 100 bias - 100 dark (600 ADU over 60s, scaled to 10s) = 900, normalised to 1.
    assert np.allclose(m_flat.data, 1.0)
    assert m_flat.meta["NCOMBINE"] == 3


def test_flats_are_normalised_before_combining(tmp_path):
    """Twilight flats fade as they are taken; the combine must see response, not level."""
    cfg = load_config()
    instrument = FixtureProfile()

    # Same 2x response pattern in every frame, at wildly different light levels.
    paths = []
    for i, level in enumerate((10000.0, 20000.0, 40000.0)):
        header = fits.Header()
        header.update({"IMAGETYP": "Flat", "EXPTIME": 3.0, "EGAIN": 1.0,
                       "READNOIS": 5.0, "FILTER": "R", "TELESCOP": "T32"})
        data = np.full(SHAPE, level, dtype=np.float32)
        data[:, :8] = level * 0.5          # half-sensitivity half of the detector
        path = tmp_path / f"flat{i}.fits"
        fits.PrimaryHDU(data=data, header=header).writeto(str(path))
        paths.append(str(path))

    m_flat = build_master_flat(paths, instrument, None, None, cfg, LOGGER)

    # The response ratio must survive; combining raw ADU would let the 40000 ADU
    # frame dominate and would put the master nowhere near 1.0.
    bright, dim = np.median(m_flat.data[:, 8:]), np.median(m_flat.data[:, :8])
    assert np.isclose(bright / dim, 2.0, rtol=1e-3)
    assert np.isclose(np.median(m_flat.data), 1.0, rtol=1e-6)

    # And the uncertainty must describe pixel noise, not the 4x level spread.
    assert np.median(m_flat.uncertainty.array) < 0.05, "stack scatter tracked the level drift"


# --- Bad-pixel mask -----------------------------------------------------------

def _flats(value=1.0):
    return {"R": _master(np.full(SHAPE, value))}


def test_bpm_flags_dead_pixels_from_flats():
    flat = np.ones(SHAPE)
    flat[2, 2] = 0.1   # dead
    flat[3, 3] = 2.0   # stuck high
    cfg = load_config()

    bpm = _build_bpm({"R": _master(flat)}, None, None, cfg, LOGGER)

    assert bpm[2, 2] and bpm[3, 3]
    assert bpm.sum() == 2


def test_bpm_flags_hot_pixels_the_flats_cannot_see():
    """Hot pixels are a dark-current defect and are invisible in a short flat."""
    dark = np.full(SHAPE, 10.0)   # 0.1 e-/s over 100s at gain 1
    dark[5, 5] = 1000.0           # 10 e-/s -- 100x the median
    cfg = load_config()
    m_dark = _master(dark, exptime=100.0, gain=1.0, bias_subtracted=True)

    bpm = _build_bpm(_flats(), m_dark, None, cfg, LOGGER)

    assert bpm[5, 5]
    assert bpm.sum() == 1, "a flat-only mask would have found nothing here"


def test_bpm_flags_unstable_pixels_from_the_bias_stack():
    scatter = np.full(SHAPE, 5.0)
    scatter[7, 7] = 500.0
    cfg = load_config()
    m_bias = _master(np.zeros(SHAPE), ncombine=10, scatter=scatter)

    bpm = _build_bpm(_flats(), None, m_bias, cfg, LOGGER)

    assert bpm[7, 7]
    assert bpm.sum() == 1


def test_bpm_skips_tests_whose_master_is_missing():
    """No darks or biases must degrade to the flat-only mask, not crash."""
    cfg = load_config()
    bpm = _build_bpm(_flats(), None, None, cfg, LOGGER)
    assert bpm is not None and bpm.sum() == 0


def test_bpm_factor_of_zero_disables_a_test():
    dark = np.full(SHAPE, 10.0)
    dark[5, 5] = 1000.0
    cfg = load_config()
    cfg.phase1.bpm_dark_rate_factor = 0.0

    bpm = _build_bpm(_flats(), _master(dark, exptime=100.0, gain=1.0), None, cfg, LOGGER)

    assert bpm.sum() == 0


def test_stack_scatter_needs_enough_frames():
    assert _stack_scatter(_master(np.zeros(SHAPE), ncombine=2, scatter=np.ones(SHAPE))) is None
    assert _stack_scatter(_master(np.zeros(SHAPE), ncombine=5, scatter=np.ones(SHAPE))) is not None


def test_outlier_threshold_survives_a_near_zero_median():
    """A very clean detector must not have its whole frame flagged."""
    rng = np.random.default_rng(0)
    values = rng.normal(1e-4, 1e-3, 10000)  # median ~0, real scatter

    threshold = _outlier_threshold(values, factor=20.0, n_sigma=5.0)

    flagged = np.count_nonzero(values > threshold)
    assert flagged < 10, f"sigma floor failed: {flagged} of 10000 flagged"


# --- Instrument-profile driven saturation -------------------------------------

def test_profile_saturation_takes_precedence_over_config(tmp_path):
    """A detector that knows its own full well should not need a config file."""
    cfg = load_config()
    cfg.phase1.saturation_adu = 50000.0
    path = _write(tmp_path, "sci.fits", 30000.0, "Light", 60.0)

    class KnownDetector(FixtureProfile):
        def get_saturation(self, header):
            return 20000.0

    assert not load_standardized_ccds(path, FixtureProfile(), cfg)[0].sat_mask.any()
    assert load_standardized_ccds(path, KnownDetector(), cfg)[0].sat_mask.all()


def test_saturation_read_from_the_header(tmp_path):
    cfg = load_config()
    path = _write(tmp_path, "sci.fits", 30000.0, "Light", 60.0)
    with fits.open(path, mode="update") as hdul:
        hdul[0].header["SATURATE"] = 25000.0

    std_ccd = load_standardized_ccds(path, FixtureProfile(), cfg)[0]

    assert std_ccd.sat_mask.all()


# --- Subframe / ROI geometry --------------------------------------------------

def test_crop_master_returns_it_unchanged_when_shapes_agree():
    master = _master(np.zeros(SHAPE))
    assert crop_master_to_frame(master, SHAPE, (0, 0), "bias", LOGGER) is master


def test_crop_master_trims_a_full_frame_master_to_the_roi():
    """A windowed camera is a routine setup, not a reason to drop the frame."""
    data = np.arange(16 * 16, dtype=float).reshape(16, 16)
    master = _master(data)
    master.uncertainty = StdDevUncertainty(np.full(SHAPE, 3.0))

    cropped = crop_master_to_frame(master, (4, 6), (3, 2), "bias", LOGGER)

    assert cropped.data.shape == (4, 6)
    assert np.array_equal(cropped.data, data[2:6, 3:9])
    # The error plane must come along, or the science ERR plane is built on sand.
    assert cropped.uncertainty.array.shape == (4, 6)


def test_crop_master_refuses_an_irreconcilable_geometry():
    """Different binning cannot be cropped away, and must not be faked."""
    master = _master(np.zeros((8, 8)))
    assert crop_master_to_frame(master, (16, 16), (0, 0), "bias", LOGGER) is None


def test_crop_master_handles_a_plain_bad_pixel_mask():
    bpm = np.zeros((16, 16), dtype=bool)
    bpm[5, 7] = True
    cropped = crop_master_to_frame(bpm, (4, 4), (6, 4), "bad-pixel mask", LOGGER)
    assert cropped.shape == (4, 4)
    assert cropped[1, 1]  # the flagged pixel, in ROI coordinates


def test_subframed_science_frame_is_reduced_against_full_frame_masters(tmp_path):
    """End to end: the ROI frame survives ISR instead of being skipped."""
    cfg = load_config()
    path = tmp_path / "sub.fits"
    header = fits.Header({"IMAGETYP": "Light", "EXPTIME": 60.0, "EGAIN": 1.0,
                          "READNOIS": 5.0, "FILTER": "R",
                          "XORGSUBF": 4, "YORGSUBF": 2})
    fits.PrimaryHDU(np.full((8, 8), 1500.0, dtype=np.float32),
                    header=header).writeto(str(path))

    std_ccds = load_standardized_ccds(str(path), FixtureProfile(), cfg)
    assert std_ccds[0].meta["subframe_origin"] == (4, 2)

    processor = UniversalProcessor(
        master_bias=_master(np.full(SHAPE, 500.0)),
        master_dark=_master(np.full(SHAPE, 60.0), exptime=60.0, bias_subtracted=True),
        master_flats={"R": _master(np.ones(SHAPE))},
        bpm=np.zeros(SHAPE, dtype=bool),
        config=cfg, logger=LOGGER,
    )
    frames = processor.process_science_frame(std_ccds)

    assert len(frames) == 1
    assert frames[0].ccd.data.shape == (8, 8)
    assert frames[0].dq.shape == (8, 8)
    # 1500 - 500 bias - 60 dark, flat of 1.0, gain 1.0 e-/ADU
    assert np.allclose(frames[0].ccd.data, 940.0)


def test_geometry_mismatch_without_an_roi_is_skipped_not_crashed(tmp_path):
    cfg = load_config()
    path = _write(tmp_path, "big.fits", 1500.0, "Light", 60.0)  # 16x16
    std_ccds = load_standardized_ccds(path, FixtureProfile(), cfg)

    processor = UniversalProcessor(
        master_bias=_master(np.full((8, 8), 500.0)), config=cfg, logger=LOGGER,
    )
    assert processor.process_science_frame(std_ccds) == []


# --- Prior-calibration guard --------------------------------------------------

def _run_dir(tmp_path, **science_cards):
    raw = tmp_path / "raw"
    raw.mkdir()
    for i in range(2):
        _write(raw, f"bias_{i}.fits", 500.0, "Bias", 0.0)
    for i in range(2):
        header = fits.Header({"IMAGETYP": "Light", "EXPTIME": 60.0, "EGAIN": 1.0,
                              "READNOIS": 5.0, "FILTER": "R"})
        header.update(science_cards)
        fits.PrimaryHDU(np.full(SHAPE, 1500.0, dtype=np.float32),
                        header=header).writeto(str(raw / f"sci_{i}.fits"))
    return str(raw), str(tmp_path / "out")


def test_precalibrated_science_frames_are_skipped(tmp_path):
    """Reducing a reduced frame is silently wrong, so it is refused by default."""
    raw, out = _run_dir(tmp_path, CALSTAT="BDF")
    run_phase1(raw, out, config=load_config(), logger=LOGGER)
    assert glob.glob(os.path.join(out, "calibrated_*.fits")) == []


def test_precalibrated_frames_can_be_opted_into(tmp_path):
    raw, out = _run_dir(tmp_path, CALSTAT="BDF")
    cfg = load_config()
    cfg.phase1.allow_precalibrated = True
    run_phase1(raw, out, config=cfg, logger=LOGGER)
    assert len(glob.glob(os.path.join(out, "calibrated_*.fits"))) == 2


def test_raw_science_frames_are_reduced_normally(tmp_path):
    """The guard must not fire on ordinary raw frames."""
    raw, out = _run_dir(tmp_path)
    run_phase1(raw, out, config=load_config(), logger=LOGGER)
    assert len(glob.glob(os.path.join(out, "calibrated_*.fits"))) == 2


# --- Scientific hardening (WP3) ----------------------------------------------

def test_the_error_budget_excludes_the_bias_pedestal():
    """Poisson noise comes from collected charge, not from an electronic offset.

    Counting a 500 ADU pedestal as signal injects tens of electrons of
    fictitious noise into every pixel, swamping the read noise on a sky-limited
    frame.
    """
    from cassa_photometry.phase1_calibration.data_models import _poisson_plus_read_noise

    gain, read_noise, pedestal = 1.4, 7.0, 500.0
    data = np.full((4, 4), pedestal + 100.0)

    with_pedestal = _poisson_plus_read_noise(data, gain, read_noise, bias_level=0.0)
    without = _poisson_plus_read_noise(data, gain, read_noise, bias_level=pedestal)

    assert without.mean() < with_pedestal.mean()
    expected = np.sqrt(100.0 * gain + read_noise**2) / gain
    assert without.mean() == pytest.approx(expected, rel=1e-6)


def test_a_negative_measurement_still_carries_the_read_noise():
    from cassa_photometry.phase1_calibration.data_models import _poisson_plus_read_noise

    sigma = _poisson_plus_read_noise(np.full((2, 2), 400.0), 1.4, 7.0, bias_level=500.0)
    assert np.all(sigma == pytest.approx(7.0 / 1.4))


def test_the_median_uncertainty_is_wider_than_the_means():
    """mad_std/sqrt(N) is the SEM of the mean; a median's is wider by sqrt(pi/2),
    and understating it understates every science frame's ERR."""
    from cassa_photometry.phase1_calibration.pipeline import MEDIAN_SEM_FACTOR

    assert MEDIAN_SEM_FACTOR == pytest.approx(1.2533, abs=1e-4)


def test_a_tiny_stack_does_not_claim_perfect_knowledge(tmp_path):
    """mad_std over one frame is exactly zero, which would make the master
    contribute nothing at all to the error budget."""
    from cassa_photometry.phase1_calibration.pipeline import _combine_with_uncertainty

    rng = np.random.default_rng(0)
    ccds = []
    for _ in range(2):
        ccd = CCDData(rng.normal(500.0, 10.0, (8, 8)), unit=u.adu)
        ccd.uncertainty = StdDevUncertainty(np.full((8, 8), 10.0))
        ccds.append(ccd)

    master = _combine_with_uncertainty(ccds)
    assert np.all(master.uncertainty.array > 0)


def test_calibration_stacks_reject_outliers():
    from cassa_photometry.phase1_calibration.pipeline import _combine_with_uncertainty

    rng = np.random.default_rng(1)
    ccds = [CCDData(rng.normal(500.0, 3.0, (16, 16)), unit=u.adu) for _ in range(6)]
    ccds[0].data[8, 8] = 50000.0          # a cosmic ray in a calibration frame

    master = _combine_with_uncertainty(ccds, sigma_clip=3.0)
    assert master.meta["NREJECT"] > 0
    assert master.data[8, 8] < 600.0


def test_vignetting_is_not_a_bad_pixel():
    """A flat threshold against the global median condemns the corners of any
    strongly vignetted optical train."""
    from cassa_photometry.phase1_calibration.pipeline import _flat_relative_response

    y, x = np.mgrid[0:128, 0:128]
    radius = np.hypot(y - 63.5, x - 63.5) / np.hypot(63.5, 63.5)
    # 0.65 vignetting puts the corners near 0.35, below the 0.5 bad-pixel cut,
    # while they are perfectly good pixels behind a heavily vignetted train.
    flat = 1.0 - 0.65 * radius**2
    flat[64, 64] = 0.05                     # one genuinely dead pixel

    assert (flat < 0.5).sum() > 100, "the test flat should look bad to a global cut"
    relative = _flat_relative_response(flat, smooth_px=31)
    bad = relative < 0.5
    assert bad.sum() < 20, "vignetting still being flagged as defects"
    assert bad[64, 64], "the genuinely dead pixel was missed"


# --- Cosmic-ray mask vetting --------------------------------------------------
#
# astroscrappy cannot tell a stellar peak from a cosmic ray in a crowded field,
# and it *replaces* what it flags. The pipeline cannot fix the detection, but it
# can recognise an answer no cosmic-ray rate could produce.

def test_a_plausible_cosmic_ray_mask_is_kept():
    from cassa_photometry.phase1_calibration.processor import implausible_cosmic_rays

    mask = np.zeros((100, 100), dtype=bool)
    mask[:5, 0] = True                      # 0.05% of the frame: an ordinary night
    assert implausible_cosmic_rays(mask, 0.01) == 0.0


def test_a_mask_covering_the_field_is_rejected():
    from cassa_photometry.phase1_calibration.processor import implausible_cosmic_rays

    mask = np.zeros((100, 100), dtype=bool)
    mask[:12] = True                        # 12%, as a globular cluster produces
    assert implausible_cosmic_rays(mask, 0.01) == pytest.approx(0.12)


def test_the_veto_can_be_switched_off():
    from cassa_photometry.phase1_calibration.processor import implausible_cosmic_rays

    mask = np.ones((10, 10), dtype=bool)
    assert implausible_cosmic_rays(mask, 1.0) == 0.0


def test_an_implausible_mask_leaves_the_science_data_untouched(tmp_path):
    """The repair is what destroys data, so a vetoed mask must not be applied."""
    cfg = load_config()
    cfg.phase1.cr_max_fraction = 0.01

    # A field of sharp peaks: astroscrappy flags a large share of it.
    rng = np.random.default_rng(0)
    data = rng.poisson(400.0, (64, 64)).astype(np.float32)
    data[::2, ::2] += 4000.0
    path = tmp_path / "crowded.fits"
    header = fits.Header({"IMAGETYP": "Light", "EXPTIME": 60.0, "EGAIN": 1.0,
                          "READNOIS": 5.0, "FILTER": "R"})
    fits.PrimaryHDU(data, header=header).writeto(str(path))

    std_ccds = load_standardized_ccds(str(path), FixtureProfile(), cfg)
    before = std_ccds[0].ccd.data.copy()
    frames = UniversalProcessor(config=cfg, logger=LOGGER).process_science_frame(std_ccds)

    frame = frames[0]
    gain = FixtureProfile().get_gain(header) or cfg.phase1.fallback_gain
    assert frame.ccd.header.get("CRVETO"), "the veto did not fire on a peaky field"
    assert not (frame.dq & DQ_COSMIC_RAY).any(), "vetoed flags must not survive"
    np.testing.assert_allclose(frame.ccd.data, before * gain, rtol=1e-6)


def test_a_normal_frame_keeps_its_cosmic_ray_repair(tmp_path):
    """The guard must not disarm cosmic-ray rejection on ordinary data."""
    cfg = load_config()
    rng = np.random.default_rng(1)
    data = rng.normal(400.0, 3.0, (64, 64)).astype(np.float32)
    data[30, 30] = 9000.0                   # one sharp hit, as a cosmic ray is
    path = tmp_path / "ordinary.fits"
    header = fits.Header({"IMAGETYP": "Light", "EXPTIME": 60.0, "EGAIN": 1.0,
                          "READNOIS": 5.0, "FILTER": "R"})
    fits.PrimaryHDU(data, header=header).writeto(str(path))

    std_ccds = load_standardized_ccds(str(path), FixtureProfile(), cfg)
    frames = UniversalProcessor(config=cfg, logger=LOGGER).process_science_frame(std_ccds)

    frame = frames[0]
    assert not frame.ccd.header.get("CRVETO"), "a single hit is not implausible"
    assert (frame.dq & DQ_COSMIC_RAY).any(), "the real cosmic ray went unflagged"
    assert frame.ccd.data[30, 30] < 9000.0, "the cosmic ray was not repaired"
