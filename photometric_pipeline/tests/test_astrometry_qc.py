"""The astrometric residual, measured the same way on every backend.

``ASTRMS`` is not decoration: phase 3 sizes its cross-match radius from it, so a
frame without one falls back to a fixed 2" and matches more loosely. The
backends disagree about what they return -- ``solve-field`` writes a matched
star table, ASTAP writes none -- so if the residual came from the solver, the
photometry would depend on which solver happened to be installed.

These tests pin the property that makes the two install paths equivalent.
"""

import numpy as np
import pytest

from cassa_photometry.phase2_integration.astrometry_qc import pair_by_position
from cassa_photometry.phase2_integration.solvers.base import SolveResult


def test_identical_positions_have_no_residual():
    ra = np.array([10.0, 10.1, 10.2])
    dec = np.array([20.0, 20.1, 20.2])
    matched = pair_by_position(ra, dec, ra, dec)
    rms_ra, rms_dec, n = SolveResult(True, matched=matched).residuals()
    assert n == 3
    assert rms_ra == pytest.approx(0.0, abs=1e-6)
    assert rms_dec == pytest.approx(0.0, abs=1e-6)


def test_a_known_offset_is_recovered_in_arcseconds():
    """One arcsecond of declination offset must read as 1", not 1 degree."""
    ra = np.array([10.0, 10.05, 10.1, 10.15])
    dec = np.array([20.0, 20.05, 20.1, 20.15])
    offset = 1.0 / 3600.0
    matched = pair_by_position(ra, dec + offset, ra, dec)
    _, rms_dec, n = SolveResult(True, matched=matched).residuals()
    assert n == 4
    assert rms_dec == pytest.approx(1.0, abs=0.02)


def test_right_ascension_is_deprojected_by_declination():
    """An RA offset spans less sky away from the equator. Without the cos(dec)
    term the residual is overstated -- by a factor of two at +60 degrees."""
    dec_value = 60.0
    ra = np.array([10.0, 10.05, 10.1, 10.15])
    dec = np.full(ra.shape, dec_value)
    # One arcsecond ON THE SKY needs a larger RA step at this declination.
    d_ra_deg = (1.0 / 3600.0) / np.cos(np.radians(dec_value))
    matched = pair_by_position(ra + d_ra_deg, dec, ra, dec)
    rms_ra, _, _ = SolveResult(True, matched=matched).residuals()
    assert rms_ra == pytest.approx(1.0, abs=0.02)


def test_pairing_is_mutual():
    """Two field sources near one catalog star must not both claim it; that is
    how a crowded field manufactures a falsely small residual."""
    field_ra = np.array([10.0, 10.00005, 50.0, 50.0001])
    field_dec = np.array([20.0, 20.0, 30.0, 30.0])
    ref_ra = np.array([10.0, 50.0])
    ref_dec = np.array([20.0, 30.0])

    matched = pair_by_position(field_ra, field_dec, ref_ra, ref_dec,
                               radius_arcsec=5.0)
    assert matched is not None
    # Each catalog star is used at most once.
    assert len(set(zip(matched[2], matched[3], strict=True))) == len(matched[2])


def test_sources_beyond_the_radius_do_not_match():
    field_ra = np.array([10.0, 11.0])
    field_dec = np.array([20.0, 21.0])
    ref_ra = np.array([10.0 + 10.0 / 3600.0, 11.0 + 10.0 / 3600.0])
    ref_dec = np.array([20.0, 21.0])
    assert pair_by_position(field_ra, field_dec, ref_ra, ref_dec,
                            radius_arcsec=2.0) is None


def test_too_few_matches_is_none_rather_than_a_fake_number():
    """One matched star gives an RMS of exactly zero, which would look like a
    perfect solve. Refusing to report is the honest answer."""
    assert pair_by_position(np.array([10.0]), np.array([20.0]),
                            np.array([10.0]), np.array([20.0])) is None


def test_empty_inputs_are_handled():
    assert pair_by_position([], [], [1.0], [2.0]) is None
    assert pair_by_position([1.0], [2.0], [], []) is None


def test_measure_returns_none_without_a_celestial_wcs(tmp_path):
    """A frame with no WCS must degrade to "no residual", never raise."""
    from astropy.io import fits

    from cassa_photometry.config import load_config
    from cassa_photometry.phase2_integration import astrometry_qc

    frame = tmp_path / "m.fits"
    fits.PrimaryHDU(np.zeros((8, 8), dtype=np.float32)).writeto(frame)
    header = fits.getheader(frame)
    assert astrometry_qc.measure(str(frame), header, load_config()) is None
