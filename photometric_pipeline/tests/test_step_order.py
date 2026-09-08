"""Reordering: what the plan resolver permits, refuses, and why.

The design rule is in ``cassa_photometry.steps``: a step that requires a token
must follow whichever included step provides it. Everything here is a
consequence of that, and the refusals matter as much as the permissions -- a
resolver that accepted any order would be worse than none, because it would let
a user produce a plausible, silently wrong reduction and believe it was checked.
"""

import pytest

from cassa_photometry.config import load_config
from cassa_photometry.steps import (
    REGISTRY,
    StepPlanError,
    describe_plan,
    resolve_plan,
    resolved_names,
)

PHASE1_DEFAULT = ["linearity", "overscan", "bad_pixel_mask", "bias", "dark",
                  "flat", "cosmic_rays", "measure_fwhm"]


def _steps(phase, **kwargs):
    config = load_config()
    toggles = getattr(config, phase).steps
    for key, value in kwargs.items():
        setattr(toggles, key, value)
    return toggles


# --- The default is unchanged -------------------------------------------------

def test_the_default_plan_is_the_canonical_order():
    """A user who writes no config must see exactly today's behaviour."""
    assert resolved_names("phase1", _steps("phase1")) == PHASE1_DEFAULT


@pytest.mark.parametrize("phase", sorted(REGISTRY))
def test_every_default_plan_satisfies_its_own_dependencies(phase):
    """The shipped order must itself be a valid topological sort -- otherwise
    the graph is describing something other than the pipeline."""
    assert resolved_names(phase, _steps(phase))


# --- Reorderings that are genuine choices -------------------------------------

def test_cosmic_rays_may_run_before_flat_fielding():
    """Both orders are defensible, so no token links them. Over-constraining
    would refuse a sound reduction."""
    order = ["linearity", "overscan", "bad_pixel_mask", "bias", "dark",
             "cosmic_rays", "flat", "measure_fwhm"]
    assert resolved_names("phase1", _steps("phase1", order=order)) == order


def test_background_subtraction_may_follow_registration():
    order = ["align", "subtract_background", "stack", "solve_wcs",
             "measure_fwhm", "visual_qa"]
    assert resolved_names("phase2", _steps("phase2", order=order)) == order


def test_phase4_stages_are_freely_reorderable():
    """They inspect four different products and share no state."""
    order = ["stage_3_photometry", "stage_0_raw", "stage_2_master",
             "stage_1_calibrated"]
    assert resolved_names("phase4", _steps("phase4", order=order)) == order


# --- Reorderings that are refused ---------------------------------------------

def test_flat_before_bias_is_refused():
    """The headline case: it produces a plausible, silently wrong image."""
    order = ["linearity", "overscan", "bad_pixel_mask", "flat", "bias", "dark",
             "cosmic_rays", "measure_fwhm"]
    with pytest.raises(StepPlanError, match="needs what that step produces"):
        resolve_plan("phase1", _steps("phase1", order=order))


def test_stacking_before_registration_is_refused():
    order = ["subtract_background", "stack", "align", "solve_wcs",
             "measure_fwhm", "visual_qa"]
    with pytest.raises(StepPlanError, match="'stack'.*before.*'align'"):
        resolve_plan("phase2", _steps("phase2", order=order))


def test_classification_before_psf_photometry_is_refused():
    """CLASS is defined as MAG_PSF minus MAG_AUTO, so it cannot precede the
    step that produces MAG_PSF."""
    order = ["aperture_correction", "zero_point", "flux_calibration", "catalog",
             "classification", "psf_photometry"]
    with pytest.raises(StepPlanError):
        resolve_plan("phase3", _steps("phase3", order=order))


def test_a_nested_step_cannot_be_resequenced():
    """PSF photometry runs *inside* catalog generation. It can be switched off,
    but there is no seam to move it to, and silently ignoring an order is the
    exact defect this design removes."""
    order = ["aperture_correction", "psf_photometry", "zero_point",
             "flux_calibration", "catalog", "classification"]
    with pytest.raises(StepPlanError, match="runs inside another step"):
        resolve_plan("phase3", _steps("phase3", order=order))


# --- Excluding is always allowed ----------------------------------------------

def test_excluding_a_step_never_violates_a_dependency():
    """Frames flat-fielded elsewhere are a real workflow. If nothing provides a
    token, nothing requires it."""
    names = resolved_names("phase1", _steps("phase1", exclude=["bias"]))
    assert "bias" not in names
    assert "dark" in names and "flat" in names


def test_exclude_and_the_boolean_form_agree():
    by_list = resolved_names("phase1", _steps("phase1", exclude=["flat"]))
    by_bool = resolved_names("phase1", _steps("phase1", flat=False))
    assert by_list == by_bool


def test_skipped_reports_both_forms():
    toggles = _steps("phase1", exclude=["flat"], cosmic_rays=False)
    assert toggles.skipped() == ["cosmic_rays", "flat"]


# --- Bad input is named precisely ---------------------------------------------

def test_an_unknown_step_in_order_is_reported():
    with pytest.raises(StepPlanError, match="unknown step"):
        resolve_plan("phase1", _steps("phase1", order=["bias", "nonsense"]))


def test_an_order_missing_a_running_step_is_reported():
    with pytest.raises(StepPlanError, match="does not mention"):
        resolve_plan("phase1", _steps("phase1", order=["bias", "dark"]))


def test_an_order_naming_a_switched_off_step_is_reported():
    order = list(PHASE1_DEFAULT)
    with pytest.raises(StepPlanError, match="switched off"):
        resolve_plan("phase1", _steps("phase1", flat=False, order=order))


def test_a_repeated_step_in_order_is_reported():
    order = ["linearity", "linearity", "overscan", "bad_pixel_mask", "bias",
             "dark", "flat", "cosmic_rays", "measure_fwhm"]
    with pytest.raises(StepPlanError, match="does not mention|repeats"):
        resolve_plan("phase1", _steps("phase1", order=order))


def test_an_invalid_plan_fails_at_config_load(tmp_path):
    """A misordered plan discovered three hours into a reduction is a wasted
    night, so it is caught when the config is read."""
    path = tmp_path / "bad.yaml"
    path.write_text(
        "phase1:\n  steps:\n    order: [flat, bias, dark, linearity, "
        "overscan, bad_pixel_mask, cosmic_rays, measure_fwhm]\n")
    with pytest.raises(StepPlanError):
        load_config(str(path))


# --- Custom steps -------------------------------------------------------------

def _write_step(tmp_path, body):
    path = tmp_path / "my_step.py"
    path.write_text(body)
    return path


def test_a_custom_step_joins_the_plan(tmp_path):
    path = _write_step(tmp_path, "def extra(**kwargs):\n    return None\n")
    toggles = _steps("phase3", custom={"extra": f"{path}:extra"})
    assert "extra" in resolved_names("phase3", toggles)


def test_a_custom_step_can_declare_dependencies(tmp_path):
    path = _write_step(
        tmp_path,
        "def extra(**kwargs):\n    return None\n"
        "extra.requires = ('catalog',)\n")
    toggles = _steps(
        "phase3", custom={"extra": f"{path}:extra"},
        order=["aperture_correction", "zero_point", "extra", "flux_calibration",
               "catalog", "psf_photometry", "classification"])
    with pytest.raises(StepPlanError, match="needs what that step produces"):
        resolve_plan("phase3", toggles)


def test_a_custom_step_may_not_shadow_a_builtin(tmp_path):
    path = _write_step(tmp_path, "def catalog(**kwargs):\n    return None\n")
    toggles = _steps("phase3", custom={"catalog": f"{path}:catalog"})
    with pytest.raises(StepPlanError, match="already"):
        resolve_plan("phase3", toggles)


def test_a_missing_custom_step_file_is_reported(tmp_path):
    toggles = _steps("phase3", custom={"extra": f"{tmp_path}/nope.py:extra"})
    with pytest.raises(StepPlanError, match="no such file"):
        resolve_plan("phase3", toggles)


def test_a_custom_step_spec_must_name_a_callable(tmp_path):
    path = _write_step(tmp_path, "extra = 42\n")
    toggles = _steps("phase3", custom={"extra": f"{path}:extra"})
    with pytest.raises(StepPlanError, match="not callable"):
        resolve_plan("phase3", toggles)


def test_a_malformed_custom_step_spec_is_reported():
    toggles = _steps("phase3", custom={"extra": "no_colon_here.py"})
    with pytest.raises(StepPlanError, match="path/to/file.py:callable"):
        resolve_plan("phase3", toggles)


def test_the_shipped_example_steps_load():
    """examples/steps/my_steps.py is documentation people copy; it must work."""
    from pathlib import Path

    example = (Path(__file__).resolve().parents[1]
               / "examples" / "steps" / "my_steps.py")
    assert example.exists()
    toggles = _steps("phase3",
                     custom={"write_bright_list": f"{example}:write_bright_list"})
    assert "write_bright_list" in resolved_names("phase3", toggles)


# --- describe_plan (what --show-plan prints) ----------------------------------

def test_describe_plan_marks_what_is_off():
    text = describe_plan("phase1", _steps("phase1", flat=False))
    assert "flat" in text and "[off]" in text
    assert "1. linearity" in text


def test_describe_plan_reports_an_invalid_plan_instead_of_raising():
    """--show-plan is what a user reaches for when their config is wrong; it
    must explain the problem rather than traceback."""
    text = describe_plan("phase1", _steps("phase1", order=["bias"]))
    assert "INVALID" in text
