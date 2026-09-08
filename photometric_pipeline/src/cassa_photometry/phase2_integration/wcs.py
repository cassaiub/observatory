"""Astrometric calibration: turn a master stack into a WCS-solved frame.

The solve itself is delegated to a backend (``solve-field`` or an in-process
solver); this module decides *what to hand it* and *what to believe afterwards*.
Three of those decisions were previously wrong in ways that cost solves:

* **The frame's own pointing outranks anything else.** The old code preferred
  the previous group's solved centre over the frame's ``OBJCTRA``/``OBJCTDEC``,
  which pointed the solver at the wrong sky for every group after the first in
  any run covering more than one field.
* **A failed solve widens its hints instead of discarding them.** The old
  fallback dropped the scale and position entirely and asked the solver to
  search the sky. Widening is both faster and far likelier to succeed.
* **A header's stated scale is a hint, not a fact.** The workshop frames carry
  ``SECPIX = 0.4`` where the truth is 0.591 -- and a *wrong* scale makes a solve
  fail, where a missing one only makes it slow. So the window is generous, the
  fallback drops the scale altogether, and after every success the solved scale
  is compared against the header's claim and the disagreement reported.

``solve-field`` only understands a single-image FITS and its output drops extra
extensions, so we solve on the SCI plane and re-attach ERR/DQ afterwards. No
resampling happens, so the planes stay aligned.
"""

import os

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from cassa_photometry.astrometry_index import IndexManifest, build_store, required_indexes
from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import read_mef, write_mef
from cassa_photometry.instruments import get_profile
from cassa_photometry.phase2_integration.solvers import (
    SolveFieldSolver,
    SolveHints,
    get_solver,
    solved_pixel_scale,
)


class WCSSolver:
    """Solve a master frame's WCS, fetching the index files it needs."""

    def __init__(self, logger, config=None, instrument=None):
        self.logger = logger
        self.config = config or load_config()
        self.instrument = instrument or get_profile(
            self.config.instrument, config=self.config
        )
        self.index_dir = self.config.resolve_astrometry_index_dir()
        self.manifest = IndexManifest.default()
        self.store = build_store(self.config)
        self.solver = get_solver(logger, self.config)

        #: Centre of the last successful solve. A *last-resort* pointing hint for
        #: a frame whose own header says nothing -- never a substitute for one
        #: that does, since a run can cover several fields.
        self.global_anchor_ra = None
        self.global_anchor_dec = None

        if self.solver is not None:
            self.logger.info("Plate solver: %s", self.solver.name)
        if not len(self.manifest) and not os.path.isdir(self.index_dir):
            self.logger.warning(
                "No astrometry index manifest and no local index directory (%s). "
                "WCS solving will fail; run cassa-index-fetch or set "
                "CASSA_ASTROMETRY_INDEX.", self.index_dir,
            )

    # --- Public ---------------------------------------------------------------
    def solve(self, group, is_rescue=False, keep_temps=False):
        filepath = group.master_filepath
        if not filepath or not os.path.exists(filepath):
            return False
        if self.solver is None:
            return False

        _, err, dq, _ = read_mef(filepath)
        with fits.open(filepath) as hdul:
            header = hdul[0].header.copy()
            data = hdul[0].data

        hints = self._hints(header, data, is_rescue)
        index_paths = self._indexes_for(hints, self.config.phase2.index_scale_lo_frac)
        self.logger.info(
            "    -> Solving with %d index file(s); %s",
            len(index_paths), hints,
        )

        result = self.solver.solve(filepath, hints, index_paths)

        if not result:
            result = self._retry(filepath, hints, result)

        if not result:
            self.logger.error("    [-] WCS FAILED: %s", result.message or "no match")
            self._clean(filepath, keep_temps)
            return False

        self._finish(filepath, result, header, err, dq)
        self._clean(filepath, keep_temps)
        group.wcs_solved = True
        return True

    # --- Hints ----------------------------------------------------------------
    def _hints(self, header, data, is_rescue):
        """What we know about this frame, in the order we trust it."""
        pixel_scale = self.instrument.get_pixel_scale(header)
        if not pixel_scale:
            recorded = header.get("PIXSCALE")
            pixel_scale = float(recorded) if recorded else None

        # The frame's own pointing first. The global anchor is only a fallback,
        # and only during a rescue pass -- a run may cover several fields.
        ra_hint = _as_degrees(header.get("OBJCTRA"), header.get("RA"), is_ra=True)
        dec_hint = _as_degrees(header.get("OBJCTDEC"), header.get("DEC"), is_ra=False)
        if ra_hint is None or dec_hint is None:
            if is_rescue and self.global_anchor_ra is not None:
                ra_hint, dec_hint = self.global_anchor_ra, self.global_anchor_dec
                self.logger.info(
                    "    -> No pointing in the header; using the last solved centre."
                )

        return SolveHints(
            ra_deg=ra_hint, dec_deg=dec_hint,
            radius_deg=self.config.phase2.index_search_radius_deg,
            pixel_scale=pixel_scale,
            scale_tolerance=self.config.phase2.solve_scale_tolerance,
            naxis1=header.get("NAXIS1"), naxis2=header.get("NAXIS2"),
            noise_adu=_measure_noise(data),
            timeout_s=self.config.phase2.solve_timeout_s,
        )

    def _indexes_for(self, hints, scale_lo_frac, scale_hi_frac=None):
        """Select and fetch the index files this field needs.

        ``scale_hi_frac`` is widened by the retry, because a wrong plate scale
        corrupts the *field size* and so the upper end of the quad window too:
        the iTelescope headers claim 0.4"/px where the truth is 0.591, which
        makes a 10' field look like 6.8' and excludes precisely the indexes that
        do solve it.
        """
        entries = required_indexes(
            hints.ra_deg, hints.dec_deg, hints.pixel_scale,
            hints.naxis1, hints.naxis2, manifest=self.manifest,
            radius_deg=hints.radius_deg, scale_lo_frac=scale_lo_frac,
            scale_hi_frac=scale_hi_frac if scale_hi_frac is not None
            else self.config.phase2.index_scale_hi_frac,
            logger=self.logger,
        )
        if not entries:
            # Nothing to select from -- no manifest, or no pointing/scale. The
            # backend falls back to whatever the local directory holds, which is
            # exactly the historical behaviour.
            return []
        return self.store.ensure(entries, logger=self.logger)

    # --- Retry ----------------------------------------------------------------
    def _retry(self, filepath, hints, first):
        """Widen the hints rather than discard them, then drop the scale.

        Two escalations, because the two failure modes are different: too narrow
        a search, and a header whose stated scale is simply wrong.
        """
        phase2 = self.config.phase2
        self.logger.info("    [!] First pass failed (%s). Widening.",
                         first.message or "no match")

        wide = hints.widened(phase2.solve_scale_tolerance_wide,
                             phase2.index_blind_radius_deg)
        wide.timeout_s = phase2.solve_timeout_wide_s
        paths = self._indexes_for(
            wide, phase2.index_scale_lo_frac_wide,
            scale_hi_frac=phase2.index_scale_hi_frac
            * (1.0 + phase2.solve_scale_tolerance_wide),
        )
        result = self.solver.solve(filepath, wide, paths)
        if result:
            return result

        # A wrong scale hint is worse than none: it excludes the true scale.
        self.logger.info(
            "    [!] Still unsolved. Retrying with no scale hint -- the header's "
            "stated scale may be wrong."
        )
        scaleless = wide.without_scale()
        scaleless.timeout_s = phase2.solve_timeout_wide_s
        return self.solver.solve(filepath, scaleless, paths)

    # --- Results --------------------------------------------------------------
    def _finish(self, filepath, result, original_header, err, dq):
        """Record the solve, its residual, and any scale disagreement."""
        header = result.header
        self.global_anchor_ra = header.get("CRVAL1")
        self.global_anchor_dec = header.get("CRVAL2")

        self._record_residual(header, result)
        self._check_scale(header, original_header)
        header["WCSSOLVR"] = (result.backend, "Plate-solving backend")

        with fits.open(filepath) as hdul:
            solved_data = hdul[0].data
        write_mef(filepath, sci=solved_data, err=err, dq=dq, header=header,
                  history=f"PHASE 2: WCS solved via {result.backend}")
        self.logger.info("    [+] SUCCESS: WCS mapped (%s).", result.backend)

    def _record_residual(self, header, result):
        """Stamp the astrometric RMS, which says whether the WCS actually fits."""
        residuals = result.residuals()
        if residuals is None:
            self.logger.warning(
                "    [!] Solved, but no matched-star table: astrometric RMS not recorded."
            )
            return
        rms_ra, rms_dec, n_stars = residuals
        total = float(np.hypot(rms_ra, rms_dec))
        header["CRDER1"] = (rms_ra / 3600.0, "[deg] RMS astrometric residual, axis 1")
        header["CRDER2"] = (rms_dec / 3600.0, "[deg] RMS astrometric residual, axis 2")
        header["ASTRMS"] = (total, "[arcsec] total RMS astrometric residual")
        header["ASTNSTAR"] = (n_stars, "Stars matched against the astrometry index")
        self.logger.info(
            "    [+] Astrometric RMS: %.3f\" RA, %.3f\" Dec, %.3f\" total (%d stars)",
            rms_ra, rms_dec, total, n_stars,
        )

    def _check_scale(self, solved_header, original_header):
        """Compare the solved plate scale against what the header claimed.

        Free, and it is how a wrong ``SECPIX`` gets noticed rather than silently
        breaking every future solve on that instrument.
        """
        measured = solved_pixel_scale(solved_header)
        if not measured:
            return
        solved_header["PIXSCALE"] = (measured, "[arcsec/pixel] Solved plate scale")

        claimed = self.instrument.get_pixel_scale(original_header)
        if not claimed:
            return
        difference = abs(measured - claimed) / measured
        if difference > self.config.phase2.solve_scale_warn_frac:
            self.logger.warning(
                "    [!] The header's plate scale (%.4f\"/px) disagrees with the "
                "solved scale (%.4f\"/px) by %.0f%%. The header is wrong; a wrong "
                "scale hint makes solves FAIL. Check SECPIX/XPIXSZ/FOCALLEN, or "
                "set detector.pixel_scale_arcsec in your config.",
                claimed, measured, 100 * difference,
            )

    def _clean(self, filepath, keep_temps):
        if not keep_temps:
            SolveFieldSolver.clean_temps(filepath)


# --- Helpers ------------------------------------------------------------------

def _measure_noise(data):
    """Background noise of a frame in its own units, or None.

    Handed to the extractor instead of a hardcoded constant: the flag it feeds
    is an assumed noise level, so a fixed value is only ever right for the one
    detector it was tuned on.
    """
    if data is None:
        return None
    try:
        finite = np.asarray(data, dtype=float)
        _, _, std = sigma_clipped_stats(finite[np.isfinite(finite)], sigma=3.0)
        return float(std) if np.isfinite(std) and std > 0 else None
    except Exception:
        return None


def _as_degrees(*values, is_ra):
    """First value that parses as degrees, accepting sexagesimal.

    Acquisition software writes ``OBJCTRA`` both ways -- decimal degrees and
    ``HH MM SS`` -- and reading one as the other is a silent pointing error of
    hours.
    """
    for value in values:
        if value is None or value == "":
            continue
        if isinstance(value, (int, float)):
            return float(value)
        text = str(value).strip().replace(":", " ")
        try:
            return float(text)
        except ValueError:
            pass
        parts = text.split()
        if len(parts) >= 2:
            try:
                sign = -1.0 if parts[0].startswith("-") else 1.0
                numbers = [abs(float(p)) for p in parts[:3]]
                while len(numbers) < 3:
                    numbers.append(0.0)
                degrees = numbers[0] + numbers[1] / 60.0 + numbers[2] / 3600.0
                # An RA written sexagesimally is in hours.
                return sign * degrees * (15.0 if is_ra else 1.0)
            except ValueError:
                continue
    return None
