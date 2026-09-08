"""Aperture correction: from the flux in an aperture to a source's total flux.

An aperture never catches all of a star's light, and the fraction it misses
depends on the PSF. Two consequences drove this module:

* **The zero point and the catalog were on different flux scales.** The zero
  point was measured in a circular aperture and then applied to isophotal
  segment fluxes, which capture a different -- and *brightness-dependent* --
  fraction of a source. That is not an offset but a **tilt in the magnitude
  scale**, and there was no aperture correction anywhere in the package.
* **The correction is not a constant across the field.** An 8-inch f/5
  Newtonian has coma, so its PSF broadens off-axis; the pipeline's own measured
  FWHM is 14% larger averaged over the field than at the centre. A single scalar
  correction is right in the middle and wrong in the corners.

So the correction is measured from a curve of growth on bright, isolated,
unsaturated stars, and fitted as a low-order surface across the field when
enough stars support it. Where they do not, it falls back to a scalar **and
folds the measured spread into the zero-point uncertainty** rather than
discarding it: an honest error bar beats a precise wrong number.

For a Gaussian PSF the correction at 2xFWHM is essentially zero, which is why a
Gaussian test suite can pass while this is missing entirely. For a Moffat with
beta = 2.5 -- a realistic atmospheric PSF -- it is about 0.074 mag.
"""

import numpy as np
from astropy.stats import sigma_clipped_stats

#: Radii sampled by the curve of growth, in units of the measured FWHM.
GROWTH_RADII_FWHM = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0])

#: Fewest usable stars before a spatial fit is attempted rather than a scalar.
MIN_STARS_FOR_SURFACE = 12

#: Fewest usable stars for any measurement at all.
MIN_STARS_FOR_SCALAR = 3


class ApertureCorrection:
    """The correction from an aperture magnitude to a total magnitude.

    Positive: a magnitude measured in the aperture is *fainter* than the truth,
    so the correction is subtracted from it (or added to the zero point).
    """

    def __init__(self, value, scatter=0.0, n_stars=0, surface=None, shape=None,
                 method="none"):
        self.value = float(value)
        #: Corner-to-centre spread, folded into MAGZERR when no surface is fitted.
        self.scatter = float(scatter)
        self.n_stars = int(n_stars)
        self._surface = surface
        self._shape = shape
        self.method = method

    def at(self, x, y):
        """The correction at a position, in magnitudes."""
        if self._surface is None:
            return np.full(np.shape(x), self.value, dtype=float) if np.ndim(x) else self.value
        return _evaluate_surface(self._surface, x, y, self._shape)

    def __repr__(self):  # pragma: no cover - debugging aid
        return (f"ApertureCorrection({self.value:+.4f} mag, {self.method}, "
                f"{self.n_stars} stars)")


def measure(data, positions, fwhm, aperture_radius, shape=None,
            plateau_factor=5.0, logger=None):
    """Measure the aperture correction for a given aperture, from a curve of growth.

    Parameters
    ----------
    data : ndarray
        Background-subtracted science image.
    positions : ndarray, shape (N, 2)
        Bright, isolated, unsaturated star positions.
    fwhm : float
        Measured FWHM in pixels; the growth radii are expressed in units of it.
    aperture_radius : float
        The radius the zero point was measured in, in pixels.
    plateau_factor : float
        Radius, in FWHM, taken as "total". Beyond a few FWHM a curve of growth is
        flat and dominated by background noise, so extending it adds scatter
        rather than flux.
    """
    from photutils.aperture import CircularAperture, aperture_photometry

    positions = np.atleast_2d(np.asarray(positions, dtype=float))
    if positions.shape[0] < MIN_STARS_FOR_SCALAR or not fwhm or fwhm <= 0:
        return ApertureCorrection(0.0, method="unmeasured")

    radii = GROWTH_RADII_FWHM * float(fwhm)
    total_radius = float(plateau_factor) * float(fwhm)

    apertures = [CircularAperture(positions, r=r) for r in radii]
    apertures.append(CircularAperture(positions, r=total_radius))
    try:
        table = aperture_photometry(np.nan_to_num(data, nan=0.0), apertures)
    except Exception:
        return ApertureCorrection(0.0, method="unmeasured")

    total = np.asarray(table[f"aperture_sum_{len(radii)}"], dtype=float)
    curve = np.array([
        np.asarray(table[f"aperture_sum_{i}"], dtype=float) for i in range(len(radii))
    ])

    usable = np.isfinite(total) & (total > 0)
    if usable.sum() < MIN_STARS_FOR_SCALAR:
        return ApertureCorrection(0.0, method="unmeasured")

    fractions = curve[:, usable] / total[usable]
    return from_curve(fractions, radii, positions[usable], aperture_radius,
                      shape=shape, logger=logger)


def from_curve(fractions, radii, positions, aperture_radius, shape=None,
               logger=None):
    """Turn a measured curve of growth into an :class:`ApertureCorrection`.

    ``fractions`` is (n_radii, n_stars); ``aperture_radius`` is the radius the
    zero point was measured in.
    """
    index = int(np.argmin(np.abs(np.asarray(radii) - float(aperture_radius))))
    enclosed = np.asarray(fractions)[index]

    good = np.isfinite(enclosed) & (enclosed > 0.2) & (enclosed <= 1.05)
    if good.sum() < MIN_STARS_FOR_SCALAR:
        return ApertureCorrection(0.0, method="unmeasured")

    corrections = -2.5 * np.log10(np.clip(enclosed[good], 1e-6, None))
    mean, median, std = sigma_clipped_stats(corrections, sigma=3.0)

    star_positions = np.asarray(positions)[good]
    n_stars = int(good.sum())

    if n_stars >= MIN_STARS_FOR_SURFACE and shape is not None:
        surface = _fit_surface(star_positions, corrections, shape)
        if surface is not None:
            if logger is not None:
                logger.info(
                    "    Aperture correction: %+.4f mag at the centre, spread %.4f "
                    "across the field (surface fit, %d stars).",
                    float(_evaluate_surface(surface, shape[1] / 2, shape[0] / 2, shape)),
                    float(std), n_stars,
                )
            return ApertureCorrection(
                float(median), scatter=float(std), n_stars=n_stars,
                surface=surface, shape=shape, method="surface",
            )

    if logger is not None:
        logger.info(
            "    Aperture correction: %+.4f +/- %.4f mag (scalar, %d stars). "
            "The spread is folded into MAGZERR, because a field-dependent PSF "
            "makes a single number right at the centre and wrong in the corners.",
            float(median), float(std), n_stars,
        )
    return ApertureCorrection(float(median), scatter=float(std), n_stars=n_stars,
                              method="scalar")


def _design(x, y, shape):
    """Normalised quadratic design matrix, so the fit is conditioned."""
    ny, nx = shape
    u = (np.asarray(x, dtype=float) / max(nx, 1)) - 0.5
    v = (np.asarray(y, dtype=float) / max(ny, 1)) - 0.5
    ones = np.ones_like(u)
    return np.column_stack([ones, u, v, u * u, u * v, v * v])


def _fit_surface(positions, values, shape):
    """Least-squares quadratic surface, or None when it is not well determined."""
    try:
        design = _design(positions[:, 0], positions[:, 1], shape)
        coefficients, *_ = np.linalg.lstsq(design, np.asarray(values, dtype=float), rcond=None)
        if not np.all(np.isfinite(coefficients)):
            return None
        return coefficients
    except Exception:
        return None


def _evaluate_surface(coefficients, x, y, shape):
    scalar = np.ndim(x) == 0
    design = _design(np.atleast_1d(x), np.atleast_1d(y), shape)
    values = design @ coefficients
    return float(values[0]) if scalar else values
