"""Solver backends, and the hint handling that decides whether a solve succeeds.

Most of the value here is in the hints rather than the solving. Three of the
four solve failures found while building this were caused by what was handed to
the solver, not by the solver itself.
"""

import numpy as np
import pytest
from astropy.io import fits
from conftest import header

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.phase2_integration.solvers import (
    BACKENDS,
    InProcessSolver,
    SolveFieldSolver,
    SolveHints,
    SolveResult,
    _can_run,
    available_backends,
    get_solver,
    solved_pixel_scale,
)
from cassa_photometry.phase2_integration.solvers.astap import AstapSolver
from cassa_photometry.phase2_integration.wcs import _as_degrees, _measure_noise

# --- Hints --------------------------------------------------------------------

def test_the_scale_window_is_generous_because_headers_can_be_wrong():
    """SECPIX=0.4 against a true 0.591 is a real case; a wrong scale makes a
    solve fail, where a missing one only makes it slow."""
    hints = SolveHints(pixel_scale=0.6, scale_tolerance=0.25)
    assert hints.scale_low == pytest.approx(0.45)
    assert hints.scale_high == pytest.approx(0.75)


def test_widening_keeps_the_pointing_and_loosens_the_rest():
    hints = SolveHints(ra_deg=10.0, dec_deg=20.0, pixel_scale=0.6,
                       scale_tolerance=0.25, radius_deg=3.0)
    wide = hints.widened(1.0, 15.0)
    assert (wide.ra_deg, wide.dec_deg) == (10.0, 20.0)
    assert wide.radius_deg == 15.0
    assert wide.scale_low < hints.scale_low and wide.scale_high > hints.scale_high


def test_the_scale_can_be_dropped_entirely_as_a_last_resort():
    hints = SolveHints(ra_deg=10.0, dec_deg=20.0, pixel_scale=0.6)
    scaleless = hints.without_scale()
    assert scaleless.scale_low is None and scaleless.scale_high is None
    assert scaleless.ra_deg == 10.0, "the pointing must survive"


def test_a_scale_tolerance_never_makes_the_window_negative():
    assert SolveHints(pixel_scale=0.6, scale_tolerance=5.0).scale_low > 0


# --- Reading a pointing out of a header --------------------------------------

@pytest.mark.parametrize("value, expected", [
    (339.2669, 339.2669),
    ("339.2669", 339.2669),
    ("22 37 04.1", 339.2671),
    ("22:37:04.1", 339.2671),
])
def test_right_ascension_is_read_in_degrees_or_hours(value, expected):
    """Acquisition software writes OBJCTRA both ways, and reading one as the
    other is a pointing error of hours."""
    assert _as_degrees(value, is_ra=True) == pytest.approx(expected, abs=1e-3)


@pytest.mark.parametrize("value, expected", [
    (34.4158, 34.4158),
    ("+34 24 56.8", 34.4158),
    ("-34 24 56.8", -34.4158),
])
def test_declination_is_read_in_degrees(value, expected):
    assert _as_degrees(value, is_ra=False) == pytest.approx(expected, abs=1e-3)


def test_a_missing_pointing_is_none_not_a_guess():
    assert _as_degrees(None, "", is_ra=True) is None
    assert _as_degrees("not a coordinate", is_ra=True) is None


# --- Noise --------------------------------------------------------------------

def test_the_noise_level_is_measured_from_the_frame():
    """solve-field's --sigma is an assumed noise level in ADU, so a hardcoded
    value is only right for the detector it was tuned on. Assuming 5 where the
    truth was 56 returned 499 sources of which 28 were real."""
    rng = np.random.default_rng(0)
    assert _measure_noise(rng.normal(500.0, 12.0, (128, 128))) == pytest.approx(12.0, rel=0.1)
    assert _measure_noise(None) is None
    assert _measure_noise(np.full((8, 8), 5.0)) is None


# --- Results ------------------------------------------------------------------

def test_residuals_deproject_right_ascension():
    """An RA offset spans less sky off the equator; not de-projecting it
    overstates the residual by 1/cos(dec)."""
    dec = 60.0
    offset_deg = 1.0 / 3600.0
    result = SolveResult(True, matched=(
        np.array([10.0 + offset_deg, 10.0 + offset_deg]), np.array([dec, dec]),
        np.array([10.0, 10.0]), np.array([dec, dec]),
    ))
    rms_ra, rms_dec, n = result.residuals()
    assert rms_ra == pytest.approx(np.cos(np.radians(dec)), abs=0.01)
    assert rms_dec == pytest.approx(0.0, abs=1e-6)
    assert n == 2


def test_residuals_need_at_least_two_matched_stars():
    assert SolveResult(True, matched=None).residuals() is None
    one = (np.array([1.0]), np.array([1.0]), np.array([1.0]), np.array([1.0]))
    assert SolveResult(True, matched=one).residuals() is None


def test_a_failed_result_is_falsey():
    assert not SolveResult(False)
    assert SolveResult(True)


def test_solved_pixel_scale_accounts_for_rotation():
    """CD1_1 alone is scale*cos(theta) and reads low on a rotated field."""
    scale_deg, angle = 0.598 / 3600.0, np.radians(35.0)
    solved = header(
        CTYPE1="RA---TAN", CTYPE2="DEC--TAN", CRVAL1=10.0, CRVAL2=20.0,
        CRPIX1=512, CRPIX2=512,
        CD1_1=-scale_deg * np.cos(angle), CD1_2=scale_deg * np.sin(angle),
        CD2_1=scale_deg * np.sin(angle), CD2_2=scale_deg * np.cos(angle),
    )
    assert solved_pixel_scale(solved) == pytest.approx(0.598, abs=0.002)
    assert solved_pixel_scale(fits.Header()) is None


# --- Backend selection --------------------------------------------------------

# These two check the *machine*, not the code: they are the "can this box
# actually solve" assertion. They skip rather than fail where nothing is
# installed, so a checkout on a machine without a solver reports a clean suite
# and one honest skip -- rather than two red tests that say nothing about the
# change being made.
requires_a_solver = pytest.mark.skipif(
    not available_backends(),
    reason="no plate solver installed here; run ./install.sh to get one",
)


@requires_a_solver
def test_the_installed_backend_is_usable():
    assert available_backends()


@requires_a_solver
def test_auto_picks_something_that_can_actually_run():
    solver = get_solver(get_logger("test"), load_config())
    assert solver is not None
    assert _can_run(type(solver), load_config())


@requires_a_solver
def test_an_explicit_backend_is_honoured_when_available():
    config = load_config()
    for name in available_backends():
        config.phase2.solver = name
        assert get_solver(get_logger("test"), config).name == name


def test_an_unavailable_backend_falls_back_rather_than_failing(monkeypatch, pipeline_logs):
    """A machine without the requested backend must still solve with another.

    Every backend's availability is pinned, including ASTAP's. Leaving any of
    them to whatever the machine happens to have makes the assertion depend on
    the developer's laptop: with ASTAP installed the fallback legitimately picks
    ASTAP, which is first in the order, and the test would fail while the code
    was behaving correctly.
    """
    config = load_config()
    config.phase2.solver = "astrometry-py"
    monkeypatch.setattr(AstapSolver, "available", classmethod(lambda cls, config=None: False))
    monkeypatch.setattr(InProcessSolver, "available", classmethod(lambda cls: False))
    monkeypatch.setattr(SolveFieldSolver, "available", classmethod(lambda cls: True))

    solver = get_solver(get_logger("test"), config)
    assert solver.name == "solve-field"
    assert "not installed" in pipeline_logs.text


def test_the_fallback_order_is_astap_then_solve_field_then_in_process(monkeypatch):
    """The order is a decision, not an accident: ASTAP is the only backend that
    exists on every supported platform, so it is tried first."""
    config = load_config()
    config.phase2.solver = "auto"
    monkeypatch.setattr(AstapSolver, "available", classmethod(lambda cls, config=None: True))
    monkeypatch.setattr(SolveFieldSolver, "available", classmethod(lambda cls: True))
    monkeypatch.setattr(InProcessSolver, "available", classmethod(lambda cls: True))
    assert get_solver(get_logger("test"), config).name == "astap"

    monkeypatch.setattr(AstapSolver, "available", classmethod(lambda cls, config=None: False))
    assert get_solver(get_logger("test"), config).name == "solve-field"

    monkeypatch.setattr(SolveFieldSolver, "available", classmethod(lambda cls: False))
    assert get_solver(get_logger("test"), config).name == "astrometry-py"


def test_no_backend_at_all_says_how_to_get_one(monkeypatch, pipeline_logs):
    config = load_config()
    for backend in BACKENDS.values():
        monkeypatch.setattr(backend, "available", classmethod(lambda cls: False))
    assert get_solver(get_logger("test"), config) is None
    assert "cassa-doctor" in pipeline_logs.text


def test_an_unknown_backend_name_is_rejected():
    config = load_config()
    config.phase2.solver = "nope"
    with pytest.raises(KeyError, match="Unknown solver"):
        get_solver(get_logger("test"), config)


# --- solve-field specifics ----------------------------------------------------

def test_the_config_names_index_files_explicitly(tmp_path):
    """Explicit `index` lines are what make on-demand fetching possible, and
    stop a user with the full 5 GB set from loading all of it."""
    solver = SolveFieldSolver(get_logger("test"), load_config())
    path = solver._write_config(["/a/index-4203-14.fits", "/a/index-4204-14.fits"])
    try:
        text = open(path).read()
        assert "index /a/index-4203-14.fits" in text
        assert "autoindex" not in text
    finally:
        import os
        os.remove(path)


def test_without_a_selection_it_falls_back_to_the_configured_directory():
    solver = SolveFieldSolver(get_logger("test"), load_config())
    path = solver._write_config([])
    try:
        text = open(path).read()
        assert "add_path" in text and "autoindex" in text
    finally:
        import os
        os.remove(path)


def test_the_command_carries_the_measured_noise_and_the_scale_window():
    solver = SolveFieldSolver(get_logger("test"), load_config())
    hints = SolveHints(ra_deg=339.27, dec_deg=34.42, pixel_scale=0.598,
                       scale_tolerance=0.25, naxis1=2048, naxis2=1400,
                       noise_adu=56.1, radius_deg=3.0)
    command = solver._command("frame.fits", "out.fits", hints, "cfg")
    assert "--sigma" in command and command[command.index("--sigma") + 1].startswith("56")
    assert "--scale-low" in command and "--scale-high" in command
    assert "--ra" in command and "--radius" in command


def test_no_noise_estimate_means_no_sigma_flag():
    """Better to let the extractor estimate than to assert a wrong number."""
    solver = SolveFieldSolver(get_logger("test"), load_config())
    command = solver._command("f.fits", "o.fits", SolveHints(noise_adu=None), "cfg")
    assert "--sigma" not in command


def test_the_in_process_backend_is_not_confused_with_the_astrometry_net_bindings():
    """conda's astrometry.net package installs a module of the same name. Only
    the PyPI solver has `Solver`; mistaking one for the other calls into the
    wrong library. Neither may be installed at all, which is also fine."""
    from cassa_photometry.phase2_integration.solvers.inprocess import _astrometry_module

    try:
        import astrometry
    except ImportError:
        assert _astrometry_module() is None
        assert not InProcessSolver.available()
        return

    if hasattr(astrometry, "Solver"):
        assert _astrometry_module() is astrometry
        assert InProcessSolver.available()
    else:
        assert _astrometry_module() is None
        assert not InProcessSolver.available()


def test_the_in_process_backend_reports_why_it_cannot_run(tmp_path):
    solver = InProcessSolver(get_logger("test"), load_config())
    result = solver.solve(str(tmp_path / "nothing.fits"), SolveHints(), index_paths=[])
    assert not result
    assert result.message
