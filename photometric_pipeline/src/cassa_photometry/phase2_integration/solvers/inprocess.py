"""In-process plate solving, for machines with no ``solve-field`` binary.

This is what makes ``pip install`` sufficient. The Astrometry.net suite cannot be
installed with pip, but the PyPI ``astrometry`` package wraps the same C library
and ships binary wheels, so a plain virtual environment can still solve a WCS.

**A name collision to be careful about.** conda-forge's ``astrometry`` package --
the one that provides ``solve-field`` -- installs a Python module *also* called
``astrometry``, exposing ``astrometry.util``. It is a different piece of
software, and it has no ``Solver``. So the import here is guarded on the
attribute rather than on the module name; getting this wrong would mean calling
into whichever of the two happened to be installed.

Source extraction is ours rather than the solver's, using ``sep`` on the science
plane with the DQ mask applied. That is an improvement on shelling out: bad
pixels and cosmic rays are excluded from the star list by construction, rather
than having to be survived by the matcher.
"""

import numpy as np

from cassa_photometry.phase2_integration.solvers.base import Solver, SolveResult

#: Maximum sources handed to the matcher, brightest first.
MAX_SOURCES = 200

#: Detection threshold in sigma above the background.
DETECT_SIGMA = 5.0


def _astrometry_module():
    """The PyPI ``astrometry`` solver package, or None.

    Identified by ``Solver`` because the astrometry.net bindings share the name.
    """
    try:
        import astrometry
    except ImportError:
        return None
    return astrometry if hasattr(astrometry, "Solver") else None


class InProcessSolver(Solver):
    """Solve with the PyPI ``astrometry`` package, extracting sources ourselves."""

    name = "astrometry-py"

    @classmethod
    def available(cls):
        return _astrometry_module() is not None

    def solve(self, filepath, hints, index_paths=()):
        astrometry = _astrometry_module()
        if astrometry is None:
            return SolveResult(
                False, backend=self.name,
                message="the PyPI 'astrometry' package is not installed "
                        '(pip install "cassa-photometry[solver]")',
            )
        if not index_paths:
            return SolveResult(
                False, backend=self.name,
                message="no index files selected; this backend cannot search a directory",
            )

        stars = self._extract(filepath)
        if stars is None or len(stars) < 4:
            return SolveResult(False, backend=self.name,
                               message=f"only {0 if stars is None else len(stars)} "
                                       "sources extracted; need at least 4")

        size_hint = None
        if hints.scale_low and hints.scale_high:
            size_hint = astrometry.SizeHint(
                lower_arcsec_per_pixel=float(hints.scale_low),
                upper_arcsec_per_pixel=float(hints.scale_high),
            )
        position_hint = None
        if hints.ra_deg is not None and hints.dec_deg is not None:
            position_hint = astrometry.PositionHint(
                ra_deg=float(hints.ra_deg), dec_deg=float(hints.dec_deg),
                radius_deg=float(hints.radius_deg),
            )

        try:
            with astrometry.Solver([str(p) for p in index_paths]) as solver:
                solution = solver.solve(
                    stars=[(float(x), float(y)) for x, y in stars],
                    size_hint=size_hint,
                    position_hint=position_hint,
                    solution_parameters=astrometry.SolutionParameters(),
                )
        except Exception as exc:
            return SolveResult(False, backend=self.name, message=f"{type(exc).__name__}: {exc}")

        if not solution.has_match():
            return SolveResult(False, backend=self.name, message="no match found")

        match = solution.best_match()
        return SolveResult(True, header=self._header(match, filepath),
                           matched=self._matched(match, stars), backend=self.name)

    # -- helpers ---------------------------------------------------------------
    def _extract(self, filepath):
        """Bright sources from the science plane, with DQ-flagged pixels masked."""
        try:
            import sep
        except ImportError:
            self.logger.warning("sep is not installed; cannot extract sources in-process.")
            return None

        from cassa_photometry.fits_utils import read_mef

        sci, _, dq, _ = read_mef(filepath)
        data = np.ascontiguousarray(np.nan_to_num(sci, nan=0.0), dtype=np.float32)
        mask = None
        if dq is not None:
            mask = np.ascontiguousarray(np.asarray(dq) != 0)

        try:
            background = sep.Background(data, mask=mask)
            subtracted = data - background.back()
            sources = sep.extract(subtracted, DETECT_SIGMA, err=background.globalrms,
                                  mask=mask)
        except Exception as exc:
            self.logger.warning("Source extraction failed: %s", exc)
            return None

        if sources is None or len(sources) == 0:
            return None
        order = np.argsort(sources["flux"])[::-1][:MAX_SOURCES]
        return np.column_stack([sources["x"][order], sources["y"][order]])

    def _header(self, match, filepath):
        """A FITS WCS header from the match's own WCS fields."""
        from astropy.io import fits

        header = fits.Header()
        wcs_fields = match.wcs_fields
        for key, value in wcs_fields.items():
            # The package hands back (value, comment) pairs.
            header[key] = value if not isinstance(value, tuple) else value[0]
        return header

    def _matched(self, match, stars):
        """Field/index positions of the matched stars, for the residual."""
        try:
            index_ra = np.array([s.ra_deg for s in match.stars], dtype=float)
            index_dec = np.array([s.dec_deg for s in match.stars], dtype=float)
        except Exception:
            return None
        try:
            from astropy.io import fits
            from astropy.wcs import WCS

            header = fits.Header()
            for key, value in match.wcs_fields.items():
                header[key] = value if not isinstance(value, tuple) else value[0]
            wcs = WCS(header)
            pixels = np.array([(s.metadata.get("x", np.nan), s.metadata.get("y", np.nan))
                               for s in match.stars], dtype=float)
            if not np.isfinite(pixels).all():
                return None
            field_ra, field_dec = wcs.all_pix2world(pixels[:, 0], pixels[:, 1], 0)
        except Exception:
            return None
        return (field_ra, field_dec, index_ra, index_dec)
