"""Backwards-compatible re-export of :mod:`cassa_photometry.psf`.

The PSF measurement moved to the top level when it stopped being a diagnostics
concern -- phases 1, 2 and 3 all depend on it now. Importing from here still
works.
"""

from cassa_photometry.psf import (  # noqa: F401
    background_rms,
    detect_stars,
    estimate_fwhm,
    image_stats,
)

__all__ = ["image_stats", "background_rms", "detect_stars", "estimate_fwhm"]
