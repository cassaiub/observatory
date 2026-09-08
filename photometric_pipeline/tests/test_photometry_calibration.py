"""Phase 3 flux calibration: the numbers the pipeline exists to produce.

Most of these compare against *simulated truth* rather than against the
pipeline's own output, because the defects being fixed here all produce
plausible-looking numbers. The headline one is measured end to end in
``test_aperture_magnitudes_are_unbiased_against_truth``.
"""

import numpy as np
import pytest
from conftest import header

from cassa_photometry.config import load_config
from cassa_photometry.phase3_photometry import apcor, photsys
from cassa_photometry.phase3_photometry.engine import (
    _mutual_matches,
    _single_frame_exposure,
)
from cassa_photometry.simulate.scene import _add_moffat, enclosed_fraction

# --- Photometric systems ------------------------------------------------------

def test_johnson_bands_are_not_on_the_ab_system():
    """The Jy conversion is the AB relation, but B/V/R/I are calibrated against
    Vega magnitudes. Ignoring the offset makes B fluxes ~8% wrong."""
    assert photsys.system_of("B") == "Vega"
    assert photsys.system_of("G") == "AB"
    assert photsys.ab_offset("G") == 0.0

    b_error = 10 ** (-0.4 * photsys.ab_offset("B")) - 1
    assert abs(b_error) > 0.08


def test_an_unknown_band_gets_no_invented_offset():
    assert photsys.ab_offset("ZZ") == 0.0


# --- Aperture correction ------------------------------------------------------

def _moffat_field(fwhm=4.2, n=40, shape=(600, 600), seed=0):
    rng = np.random.default_rng(seed)
    image = np.zeros(shape)
    positions = np.column_stack([
        rng.uniform(60, shape[1] - 60, n), rng.uniform(60, shape[0] - 60, n)
    ])
    for x, y in positions:
        _add_moffat(image, x, y, 200000.0, fwhm)
    return image, positions


def test_the_aperture_correction_matches_the_analytic_moffat_value():
    """For a Gaussian this correction is ~0, which is why a Gaussian test suite
    passes while it is missing. For a realistic Moffat it is ~0.074 mag."""
    fwhm = 4.2
    image, positions = _moffat_field(fwhm)
    radius = 2.0 * fwhm

    correction = apcor.measure(image, positions, fwhm, aperture_radius=radius,
                               shape=image.shape)
    expected = -2.5 * np.log10(
        enclosed_fraction(radius, fwhm) / enclosed_fraction(5 * fwhm, fwhm)
    )
    assert correction.value == pytest.approx(expected, abs=0.005)
    assert correction.n_stars > 20


def test_a_field_dependent_correction_is_fitted_as_a_surface():
    """An 8-inch Newtonian has coma, so the PSF broadens off-axis and a scalar
    correction is right at the centre and wrong in the corners."""
    shape = (1400, 2048)
    rng = np.random.default_rng(0)
    positions = np.column_stack([rng.uniform(0, 2048, 200), rng.uniform(0, 1400, 200)])
    u = positions[:, 0] / 2048 - 0.5
    v = positions[:, 1] / 1400 - 0.5
    values = 0.05 + 0.4 * (u**2 + v**2)

    surface = apcor._fit_surface(positions, values, shape)
    assert surface is not None
    assert apcor._evaluate_surface(surface, 1024, 700, shape) == pytest.approx(0.05, abs=0.01)
    assert apcor._evaluate_surface(surface, 0, 0, shape) == pytest.approx(0.25, abs=0.01)


def test_too_few_stars_gives_no_correction_rather_than_a_wrong_one():
    image, positions = _moffat_field(n=2)
    correction = apcor.measure(image, positions[:2], 4.2, aperture_radius=8.4)
    assert correction.method == "unmeasured"
    assert correction.value == 0.0


def test_a_scalar_correction_reports_its_spread_for_the_error_budget():
    """An honest error bar beats a precise wrong number."""
    correction = apcor.ApertureCorrection(0.07, scatter=0.03, n_stars=5, method="scalar")
    assert correction.scatter == 0.03
    assert correction.at(0, 0) == correction.at(1000, 1000)


# --- Exposure normalisation ---------------------------------------------------

def test_exposure_normalisation_uses_the_single_frame_exposure():
    """A master is a weighted MEAN, so its pixels are electrons per single
    frame. Dividing by TOT_EXP would be wrong by a factor of N."""
    assert _single_frame_exposure(header(EXPMEAN=118.4, EXPTIME=120, TOT_EXP=600)) == 118.4


def test_exptime_is_the_fallback_and_tot_exp_is_never_used():
    assert _single_frame_exposure(header(EXPTIME=120, TOT_EXP=600)) == 120.0
    assert _single_frame_exposure(header(TOT_EXP=600)) == 1.0


# --- Cross-matching -----------------------------------------------------------

def test_matching_requires_agreement_in_both_directions():
    """One-directional matching lets two detections claim the same catalog star
    -- a deblended pair, or a star beside a galaxy -- and both enter the ZP."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    catalog = SkyCoord([10.0] * u.deg, [20.0] * u.deg)
    # Two detections near one catalog star: only the closer may match.
    detections = SkyCoord([10.00003, 10.00008] * u.deg, [20.0, 20.0] * u.deg)

    pairs = _mutual_matches(detections, catalog, 2.0 * u.arcsec)
    assert len(pairs) == 1
    assert pairs[0][0] == 0


def test_matches_beyond_the_radius_are_rejected():
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    catalog = SkyCoord([10.0] * u.deg, [20.0] * u.deg)
    detections = SkyCoord([10.01] * u.deg, [20.0] * u.deg)
    assert _mutual_matches(detections, catalog, 2.0 * u.arcsec) == []


def test_an_empty_catalog_matches_nothing():
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    empty = SkyCoord([] * u.deg, [] * u.deg)
    some = SkyCoord([10.0] * u.deg, [20.0] * u.deg)
    assert _mutual_matches(some, empty, 2.0 * u.arcsec) == []


# --- The match radius ---------------------------------------------------------

def test_the_zero_point_uses_the_measured_astrometric_residual():
    """The engine used a hardcoded 2 arcsec while the verification tool already
    derived one from ASTRMS; on the workshop masters that was ~4 sigma."""
    from cassa_photometry.phase3_photometry.verify import match_radius_from_header

    config = load_config()
    radius, why = match_radius_from_header(header(ASTRMS=0.5, ASTNSTAR=60), config)
    assert radius == pytest.approx(1.5)
    assert "ASTRMS" in why

    fallback, why = match_radius_from_header(header(), config)
    assert fallback == config.phase3.zp_match_tol_arcsec
    assert "no ASTRMS" in why


def test_the_match_radius_is_clamped_at_both_ends():
    """A superb solve still must not match at a few tenths of an arcsecond:
    reference positions carry their own error and stars move."""
    from cassa_photometry.phase3_photometry.verify import match_radius_from_header

    config = load_config()
    tiny, _ = match_radius_from_header(header(ASTRMS=0.01), config)
    huge, _ = match_radius_from_header(header(ASTRMS=100.0), config)
    assert tiny == config.phase3.match_radius_min_arcsec
    assert huge == config.phase3.match_radius_max_arcsec


# --- The search cone ----------------------------------------------------------

def test_the_catalog_cone_covers_the_corners_of_the_field():
    """`nx/2` is the inscribed radius, not the half-diagonal: measured on a real
    master, the old expression left half the field area with no reference stars."""
    from astropy.wcs import WCS

    from cassa_photometry.phase3_photometry.engine import UniversalPhotometryEngine

    scale = 0.598 / 3600.0
    wcs = WCS(header(CTYPE1="RA---TAN", CTYPE2="DEC--TAN", CRVAL1=10.0, CRVAL2=20.0,
                     CRPIX1=1024, CRPIX2=700, CD1_1=-scale, CD1_2=0.0,
                     CD2_1=0.0, CD2_2=scale))
    engine = UniversalPhotometryEngine(config=load_config())
    radius = engine.search_radius(wcs, (1400, 2048)).to_value("deg")

    half_diagonal = 0.5 * np.hypot(2048, 1400) * scale
    assert radius == pytest.approx(half_diagonal, rel=1e-6)
    # ...and it is meaningfully larger than the old half-width expression.
    assert radius > 1.2 * (1024 * scale)
