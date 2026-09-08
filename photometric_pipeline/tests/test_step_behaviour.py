"""Switching a step off must change the product, not merely the log.

``test_step_plan.py`` proves that every declared toggle has *a* consumer;
``test_config.py`` proves the YAML parses. Neither proves the pipeline behaves
differently, which is the only thing the user actually asked for. These tests
close that gap: each one turns a step off and asserts on the result.
"""

import numpy as np
import pytest

from cassa_photometry.config import load_config
from cassa_photometry.phase2_integration.pipeline import IntegrationPipeline


def _pipeline(tmp_path, config):
    """An IntegrationPipeline instance without touching disk or hardware."""
    return IntegrationPipeline(str(tmp_path), output_dir=str(tmp_path),
                               config=config)


def _sky_plus_star(shape=(64, 64)):
    """A frame with a strong linear sky gradient and one star on it."""
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    sky = 100.0 + 0.5 * xx + 0.25 * yy
    star = 500.0 * np.exp(-(((xx - 32) ** 2 + (yy - 32) ** 2) / (2 * 2.0 ** 2)))
    return sky + star


# --- phase2.steps.subtract_background -----------------------------------------

def test_background_is_removed_by_default(tmp_path):
    config = load_config()
    pipeline = _pipeline(tmp_path, config)
    data = _sky_plus_star()

    result, background = pipeline._remove_background(data, config.phase2)

    # The sky is gone: the frame is centred near zero where it was near 100.
    assert abs(float(np.median(result))) < 0.1 * abs(float(np.median(data)))
    assert not np.allclose(result, data)
    assert np.any(np.asarray(background) != 0)


def test_background_is_kept_when_the_step_is_off(tmp_path):
    """Extended-source work needs the sky left in -- the outer disk of a
    resolved galaxy is what a background model otherwise eats."""
    config = load_config()
    config.phase2.steps.subtract_background = False
    pipeline = _pipeline(tmp_path, config)
    data = _sky_plus_star()

    result, background = pipeline._remove_background(data, config.phase2)

    assert np.allclose(result, data)
    # The level is still measured and returned, so BKGLEVEL is recorded either
    # way. It is the subtraction that is optional, not the measurement.
    assert np.any(np.asarray(background) != 0)


# --- phase2.steps.measure_fwhm ------------------------------------------------

def test_master_fwhm_is_not_measured_when_the_step_is_off(tmp_path):
    config = load_config()
    config.phase2.steps.measure_fwhm = False
    pipeline = _pipeline(tmp_path, config)

    assert pipeline._measure_master_fwhm(_sky_plus_star()) is None


def test_an_unmeasured_master_fwhm_is_absent_rather_than_the_anchors(tmp_path):
    """The master's header is copied from the anchor *frame*, which phase 1
    stamped with its own FWHMPX. If phase 2 does not measure the stack's PSF,
    that inherited card must go: phase 3 sizes every aperture from FWHMPX, so a
    stale single-frame value is read as the stack's seeing and silently
    mis-sizes the photometry."""
    from astropy.io import fits

    from cassa_photometry.fits_utils import read_mef
    from cassa_photometry.phase2_integration.io import FITSHandler

    anchor = tmp_path / "anchor.fits"
    fits.PrimaryHDU(
        np.ones((8, 8), dtype=np.float32),
        fits.Header({"FWHMPX": 4.2, "EXPTIME": 60.0}),
    ).writeto(anchor)

    out = tmp_path / "master.fits"
    FITSHandler.save_master(
        np.ones((8, 8)), np.ones((8, 8)), np.zeros((8, 8), dtype=np.int32),
        str(anchor), str(out), 3,
        {"focal_length_mm": 0, "pixel_size_um": 0, "exposure": 60.0,
         "stack_stats": {"fwhm": None}},
    )

    assert "FWHMPX" not in read_mef(str(out))[3]


# --- phase2.steps.align -------------------------------------------------------

def test_align_off_passes_the_frame_through_unregistered(tmp_path):
    config = load_config()
    config.phase2.steps.align = False
    pipeline = _pipeline(tmp_path, config)

    reference = _sky_plus_star()
    target = _sky_plus_star() * 1.3
    frame = {"fwhm": 3.0, "path": "t.fits"}
    anchor = {"fwhm": 3.0}

    data, variance, dq, scale, n_stars, noise = pipeline._register_frame(
        target, None, None, reference, anchor, frame, config.phase2)

    assert np.allclose(data, target), "the frame must not be resampled"
    assert scale == 1.0, "no star matches means no transparency normalisation"
    assert n_stars == 0
    assert dq.shape == target.shape
    assert variance.shape == target.shape
    assert noise > 0


def test_align_off_refuses_a_frame_on_a_different_grid(tmp_path):
    """Without registration the frames must already share a pixel grid. A
    mismatch is an error, not something to silently resample away."""
    config = load_config()
    config.phase2.steps.align = False
    pipeline = _pipeline(tmp_path, config)

    with pytest.raises(ValueError, match="pixel grid"):
        pipeline._register_frame(
            _sky_plus_star((32, 32)), None, None, _sky_plus_star((64, 64)),
            {"fwhm": 3.0}, {"fwhm": 3.0, "path": "t.fits"}, config.phase2)


# --- phase4.steps -------------------------------------------------------------

def test_a_phase4_stage_can_be_switched_off():
    config = load_config()
    assert config.phase4.steps.enabled("stage_0_raw")
    config.phase4.steps.stage_0_raw = False
    assert not config.phase4.steps.enabled("stage_0_raw")
    assert config.phase4.steps.skipped() == ["stage_0_raw"]


# --- Deprecated duplicate flags ----------------------------------------------

@pytest.mark.parametrize(("phase", "legacy", "step"), [
    ("phase1", "apply_linearity", "linearity"),
    ("phase2", "subtract_background", "subtract_background"),
    ("phase3", "psf_photometry", "psf_photometry"),
    ("phase3", "aperture_correction", "aperture_correction"),
])
def test_a_legacy_flag_still_switches_its_step_off(tmp_path, phase, legacy, step):
    """These names are in existing configs, the README and the notebooks. They
    keep working -- a setting that quietly stops having an effect is the same
    defect this change exists to remove."""
    path = tmp_path / "cfg.yaml"
    path.write_text(f"{phase}:\n  {legacy}: false\n")

    config = load_config(str(path))

    assert not getattr(config, phase).steps.enabled(step)


def test_the_steps_toggle_is_enough_on_its_own(tmp_path):
    """The point of collapsing the duplicates: setting the step alone works,
    without also having to find and unset the legacy flag."""
    path = tmp_path / "cfg.yaml"
    path.write_text("phase3:\n  steps:\n    psf_photometry: false\n")

    config = load_config(str(path))

    assert not config.phase3.steps.enabled("psf_photometry")
    # The legacy field is untouched and no longer gates anything.
    assert config.phase3.psf_photometry is True


# --- Custom steps actually execute --------------------------------------------

def _fake_frame(ccd):
    """A minimal StandardCCD stand-in carrying every key the ISR stages read.

    The CCD must carry a real ``fits.Header``: a bare dict stores the
    ``(value, comment)`` tuple verbatim, which is not how the pipeline's frames
    behave and would let a header assertion pass or fail for the wrong reason.
    """
    return type("Std", (), {"ccd": ccd, "sat_mask": None, "meta": {
        "filter": "V", "exposure": 1.0, "gain": 1.0, "read_noise": 10.0,
        "fringe_needed": False, "saturation_adu": 50000.0, "bias_level": 0.0,
    }})()


def test_a_custom_phase1_step_runs_on_the_frame(tmp_path, monkeypatch):
    """The plan can carry a user's own step, and phase 1 calls it with the
    frame. Resolving it is not enough -- it has to actually run."""
    from cassa_photometry.phase1_calibration.processor import UniversalProcessor

    marker = tmp_path / "ran.txt"
    step = tmp_path / "my_step.py"
    step.write_text(
        "def note(*, ccd, **_):\n"
        f"    open({str(marker)!r}, 'w').write(str(ccd.data.shape))\n"
    )

    config = load_config()
    config.phase1.steps.custom = {"note": f"{step}:note"}
    config.phase1.steps.order = [
        "linearity", "overscan", "bad_pixel_mask", "bias", "dark", "flat",
        "cosmic_rays", "note", "measure_fwhm",
    ]

    from astropy.io import fits
    from astropy.nddata import CCDData

    ccd = CCDData(np.ones((8, 8)), unit="adu", meta=fits.Header())
    std = _fake_frame(ccd)

    processor = UniversalProcessor(config=config)
    processor.process_science_frame([std])

    assert marker.exists(), "the custom step was resolved but never called"
    assert marker.read_text() == "(8, 8)"


def test_the_resolved_plan_is_recorded_in_the_frame_header(tmp_path):
    """A reordered reduction must not be indistinguishable from a default one."""
    from astropy.io import fits
    from astropy.nddata import CCDData

    from cassa_photometry.phase1_calibration.processor import UniversalProcessor

    config = load_config()
    config.phase1.steps.order = [
        "linearity", "overscan", "bad_pixel_mask", "bias", "dark",
        "cosmic_rays", "flat", "measure_fwhm",
    ]
    ccd = CCDData(np.ones((8, 8)), unit="adu", meta=fits.Header())
    std = _fake_frame(ccd)

    frames = UniversalProcessor(config=config).process_science_frame([std])

    plan = frames[0].ccd.header["CALPLAN"].split(",")
    assert plan.index("cosmic_rays") < plan.index("flat")
