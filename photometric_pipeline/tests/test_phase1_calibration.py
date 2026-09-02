import glob
import os

import numpy as np
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.instruments import ITelescopeNetworkProfile
from cassa_photometry.phase1_calibration.data_models import load_standardized_ccds
from cassa_photometry.phase1_calibration.processor import (
    UniversalProcessor, crop_master_to_frame, subtract_scaled_dark,
)
from cassa_photometry.phase1_calibration.pipeline import (
    _build_bpm, _outlier_threshold, _stack_scatter,
    build_master_bias, build_master_dark, build_master_flat,
    run as run_phase1,
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
    instrument = ITelescopeNetworkProfile()
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
    instrument = ITelescopeNetworkProfile()
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
    instrument = ITelescopeNetworkProfile()
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
    instrument = ITelescopeNetworkProfile()
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
    instrument = ITelescopeNetworkProfile()

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

    class KnownDetector(ITelescopeNetworkProfile):
        def get_saturation(self, header):
            return 20000.0

    assert not load_standardized_ccds(path, ITelescopeNetworkProfile(), cfg)[0].sat_mask.any()
    assert load_standardized_ccds(path, KnownDetector(), cfg)[0].sat_mask.all()


def test_saturation_read_from_the_header(tmp_path):
    cfg = load_config()
    path = _write(tmp_path, "sci.fits", 30000.0, "Light", 60.0)
    with fits.open(path, mode="update") as hdul:
        hdul[0].header["SATURATE"] = 25000.0

    std_ccd = load_standardized_ccds(path, ITelescopeNetworkProfile(), cfg)[0]

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

    std_ccds = load_standardized_ccds(str(path), ITelescopeNetworkProfile(), cfg)
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
    std_ccds = load_standardized_ccds(path, ITelescopeNetworkProfile(), cfg)

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
