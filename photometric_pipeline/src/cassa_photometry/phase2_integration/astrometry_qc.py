"""Astrometric quality, measured the same way whichever backend solved.

A solve that *converged* and a solve that actually *fits* are different things,
and the difference is the residual between where the WCS puts a star and where
the sky says it is. Phase 3 depends on this: it sizes its cross-match radius
from ``ASTRMS``, so a missing residual silently widens matching to a fixed 2".

The backends disagree about what they hand back:

* ``solve-field`` writes a ``.corr`` table of matched field/index pairs.
* the in-process solver returns the catalog stars it matched, but no pixels.
* **ASTAP returns neither** -- only quad counts.

So the residual cannot be left to the backend. This module computes it from
things every backend produces: the solved WCS, the sources we can detect
ourselves, and a set of reference positions. Where a backend supplies its own
catalog stars those are used (no network); otherwise the reference catalog phase
3 already fetches and caches is reused.
"""

import numpy as np

#: Detection threshold in sigma for the sources used to measure the residual.
DETECT_SIGMA = 5.0

#: Most sources to use. The residual is a statistic, not a survey; a few hundred
#: is far past the point where adding more changes the answer.
MAX_SOURCES = 300


def pair_by_position(field_ra, field_dec, ref_ra, ref_dec, radius_arcsec=5.0):
    """Pair two position lists by nearest neighbour within a radius.

    Returns ``(field_ra, field_dec, ref_ra, ref_dec)`` of the matched pairs, in
    the shape :meth:`SolveResult.residuals` expects, or None when nothing
    matches. The pairing is **mutual**: a field source and a catalog star must
    each be the other's nearest, which is what stops one bright star in a
    crowded field from claiming several catalog entries.
    """
    field_ra = np.asarray(field_ra, dtype=float)
    field_dec = np.asarray(field_dec, dtype=float)
    ref_ra = np.asarray(ref_ra, dtype=float)
    ref_dec = np.asarray(ref_dec, dtype=float)
    if field_ra.size == 0 or ref_ra.size == 0:
        return None

    try:
        import astropy.units as u
        from astropy.coordinates import SkyCoord
    except ImportError:  # pragma: no cover - astropy is a hard dependency
        return None

    field = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg)
    ref = SkyCoord(ra=ref_ra * u.deg, dec=ref_dec * u.deg)

    idx_f, d2d_f, _ = field.match_to_catalog_sky(ref)
    idx_r, _, _ = ref.match_to_catalog_sky(field)

    keep = []
    for i, (j, sep) in enumerate(zip(idx_f, d2d_f.arcsec, strict=True)):
        if sep <= radius_arcsec and idx_r[j] == i:
            keep.append((i, j))
    if len(keep) < 2:
        return None

    fi = np.array([k[0] for k in keep])
    ri = np.array([k[1] for k in keep])
    return (field_ra[fi], field_dec[fi], ref_ra[ri], ref_dec[ri])


def detect_sources(filepath, logger=None, max_sources=MAX_SOURCES):
    """Bright sources in the science plane, with DQ-flagged pixels masked.

    The same extraction the in-process backend uses, so the residual is measured
    on the same sources whichever backend solved the frame.
    """
    try:
        import sep
    except ImportError:
        if logger:
            logger.debug("sep is not installed; cannot measure an astrometric residual.")
        return None

    from cassa_photometry.fits_utils import read_mef

    try:
        sci, _, dq, _ = read_mef(filepath)
        data = np.ascontiguousarray(np.nan_to_num(sci, nan=0.0), dtype=np.float32)
        mask = np.ascontiguousarray(np.asarray(dq) != 0) if dq is not None else None
        background = sep.Background(data, mask=mask)
        sources = sep.extract(data - background.back(), DETECT_SIGMA,
                              err=background.globalrms, mask=mask)
    except Exception as exc:
        if logger:
            logger.debug("Source extraction for the residual failed: %s", exc)
        return None

    if sources is None or len(sources) == 0:
        return None
    order = np.argsort(sources["flux"])[::-1][:max_sources]
    return np.column_stack([sources["x"][order], sources["y"][order]])


def reference_positions(header, config, logger=None):
    """Reference-catalog sky positions covering this frame, or None.

    Reuses phase 3's fetcher, so the on-disk cache is shared and ``offline``
    is honoured: a field reduced once needs no network the next time, and an
    offline run that has no cache degrades to "no residual recorded" rather
    than failing the solve.
    """
    try:
        import astropy.units as u
        from astropy.wcs import WCS
        from astropy.wcs.utils import proj_plane_pixel_scales

        from cassa_photometry.phase3_photometry.catalogs import (
            fetch_reference_catalog,
        )
    except Exception:
        return None

    try:
        wcs = WCS(header)
        if not wcs.has_celestial:
            return None
        ny = int(header.get("NAXIS2") or 0)
        nx = int(header.get("NAXIS1") or 0)
        if not nx or not ny:
            return None
        centre = wcs.pixel_to_world(nx / 2, ny / 2)
        scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial)))
        radius = (0.5 * np.hypot(nx, ny) * scale) * u.deg

        catalog = fetch_reference_catalog(centre, radius, "R", logger,
                                          config=config)
    except Exception as exc:
        if logger:
            logger.debug("No reference catalog for the astrometric residual: %s", exc)
        return None

    if catalog is None or len(catalog) == 0:
        return None
    return (np.asarray(catalog["RAJ2000"], dtype=float),
            np.asarray(catalog["DEJ2000"], dtype=float))


def measure(filepath, header, config, logger=None, radius_arcsec=5.0):
    """Matched field/reference pairs for a solved frame, or None.

    Returns the 4-tuple :meth:`SolveResult.residuals` consumes, so a backend
    that cannot produce its own matched stars still gets an ``ASTRMS``.
    """
    from astropy.wcs import WCS

    pixels = detect_sources(filepath, logger)
    if pixels is None or len(pixels) < 2:
        return None

    try:
        wcs = WCS(header)
        if not wcs.has_celestial:
            return None
        field_ra, field_dec = wcs.all_pix2world(pixels[:, 0], pixels[:, 1], 0)
    except Exception:
        return None

    reference = reference_positions(header, config, logger)
    if reference is None:
        return None

    return pair_by_position(field_ra, field_dec, reference[0], reference[1],
                            radius_arcsec=radius_arcsec)
