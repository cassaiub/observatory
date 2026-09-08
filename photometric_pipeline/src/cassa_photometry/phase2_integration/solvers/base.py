"""One interface for plate solving, whichever backend does the work.

Two backends are supported and either is sufficient:

* ``solve-field``, the Astrometry.net binary. Mature, and preferred when it is
  installed.
* an in-process solver from the PyPI ``astrometry`` package, so a pip-only
  machine with no system packages can still solve a WCS.

They share this interface so the phase does not care which ran, and so the
astrometric residual is computed the same way for both -- a solve that
*converged* and a solve that actually *fits* are different things, and the
distinction has to be measured identically or it cannot be compared.
"""

import numpy as np


class SolveHints:
    """What is known about a frame before solving it.

    Every field is optional. A hint that is absent only makes the solve slower;
    a hint that is **wrong** makes it fail, which is why the scale carries a
    tolerance rather than being taken as exact.
    """

    __slots__ = ("ra_deg", "dec_deg", "radius_deg", "pixel_scale",
                 "scale_tolerance", "naxis1", "naxis2", "noise_adu", "timeout_s")

    def __init__(self, ra_deg=None, dec_deg=None, radius_deg=3.0, pixel_scale=None,
                 scale_tolerance=0.25, naxis1=None, naxis2=None, noise_adu=None,
                 timeout_s=120.0):
        self.ra_deg = ra_deg
        self.dec_deg = dec_deg
        self.radius_deg = radius_deg
        self.pixel_scale = pixel_scale
        self.scale_tolerance = scale_tolerance
        self.naxis1 = naxis1
        self.naxis2 = naxis2
        self.noise_adu = noise_adu
        self.timeout_s = timeout_s

    @property
    def scale_low(self):
        if not self.pixel_scale:
            return None
        return self.pixel_scale * (1.0 - min(self.scale_tolerance, 0.95))

    @property
    def scale_high(self):
        if not self.pixel_scale:
            return None
        return self.pixel_scale * (1.0 + self.scale_tolerance)

    def widened(self, tolerance, radius_deg):
        """A copy with a wider scale window and search radius."""
        return SolveHints(
            ra_deg=self.ra_deg, dec_deg=self.dec_deg, radius_deg=radius_deg,
            pixel_scale=self.pixel_scale, scale_tolerance=tolerance,
            naxis1=self.naxis1, naxis2=self.naxis2, noise_adu=self.noise_adu,
            timeout_s=self.timeout_s,
        )

    def without_scale(self):
        """A copy with the scale hint dropped, for when the header is wrong."""
        return SolveHints(
            ra_deg=self.ra_deg, dec_deg=self.dec_deg, radius_deg=self.radius_deg,
            pixel_scale=None, naxis1=self.naxis1, naxis2=self.naxis2,
            noise_adu=self.noise_adu, timeout_s=self.timeout_s,
        )

    def __repr__(self):  # pragma: no cover - debugging aid
        scale = f"{self.pixel_scale:.3f}" if self.pixel_scale else "?"
        return (f"SolveHints(ra={self.ra_deg}, dec={self.dec_deg}, "
                f"scale={scale}+/-{self.scale_tolerance:.0%})")


class SolveResult:
    """The outcome of one solve attempt."""

    def __init__(self, solved, header=None, matched=None, backend=None, message=None):
        self.solved = bool(solved)
        #: The solved WCS header, or None.
        self.header = header
        #: (field_ra, field_dec, index_ra, index_dec) arrays of matched stars.
        self.matched = matched
        self.backend = backend
        self.message = message

    def __bool__(self):
        return self.solved

    def residuals(self):
        """Per-axis RMS of the star/catalog residuals in arcsec, and the count.

        Returns ``(rms_ra, rms_dec, n_stars)`` or None. This is the difference
        between a solve that converged and one that fits.
        """
        if not self.matched:
            return None
        field_ra, field_dec, index_ra, index_dec = (
            np.asarray(a, dtype=float) for a in self.matched
        )
        good = (np.isfinite(field_ra) & np.isfinite(field_dec)
                & np.isfinite(index_ra) & np.isfinite(index_dec))
        if int(good.sum()) < 2:
            return None
        # An RA offset spans less sky off the equator, so de-project it.
        cos_dec = np.cos(np.radians(index_dec[good]))
        d_ra = (field_ra[good] - index_ra[good]) * cos_dec * 3600.0
        d_dec = (field_dec[good] - index_dec[good]) * 3600.0
        return (float(np.sqrt(np.mean(d_ra ** 2))),
                float(np.sqrt(np.mean(d_dec ** 2))),
                int(good.sum()))


class Solver:
    """Base class for a plate-solving backend."""

    name = "base"

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    @classmethod
    def available(cls):
        """True when this backend can actually run on this machine."""
        return False

    def solve(self, filepath, hints, index_paths=()):  # pragma: no cover - abstract
        raise NotImplementedError


def solved_pixel_scale(header):
    """Plate scale of a solved WCS in arcsec/pixel, or None.

    Uses ``proj_plane_pixel_scales`` rather than ``CD1_1``: on a rotated field
    ``CD1_1`` is ``scale * cos(theta)`` and reads low.
    """
    try:
        from astropy.wcs import WCS
        from astropy.wcs.utils import proj_plane_pixel_scales

        wcs = WCS(header)
        if not wcs.has_celestial:
            return None
        return float(np.mean(proj_plane_pixel_scales(wcs.celestial)) * 3600.0)
    except Exception:
        return None
