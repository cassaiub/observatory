"""Plate solving with ASTAP.

ASTAP is a single self-contained binary (under 1 MB) with no Python dependency,
packaged for Debian/Ubuntu as ``astap-cli`` and published for Linux x86-64 and
aarch64 and for macOS on both Intel and Apple Silicon. That makes it the one
backend installable on every platform this pipeline supports, including ARM
Linux, where neither conda-forge nor PyPI ships an astrometry.net solver at all.

**What it does not give us.** ASTAP reports quad counts, a star count and the
pointing offset, but no matched field/catalog star pairs and -- on a sparse
field -- no SIP distortion (*"Not enough quads for calculating SIP"*). Phase 3
sizes its cross-match radius from ``ASTRMS``, so the residual is measured by
``astrometry_qc`` instead of being taken from the solver. Every backend is
treated the same way there, which is what makes the two install paths produce
the same products.

**Two failure modes, deliberately distinguished.** ASTAP answers
``ERROR=Error reading star database`` when the tiles covering the field are
absent, and ``No solution found`` when it has the sky but cannot match the
stars. The first is recoverable by fetching tiles; the second is not, and
reporting them as one error would send a user hunting for a database they
already have.
"""

import os
import shutil
import subprocess

from cassa_photometry.phase2_integration.solvers.base import Solver, SolveResult

#: Files ASTAP leaves beside its output.
TEMP_SUFFIXES = (".ini", ".wcs", ".log")

#: Candidate binary names, in the order they are preferred. Debian installs the
#: command-line build as ``astap_cli``; the upstream zip ships the same name.
BINARIES = ("astap_cli", "astap")


def find_binary(configured=None):
    """Path to an ASTAP executable, or None."""
    if configured:
        path = os.path.abspath(os.path.expanduser(configured))
        return path if os.path.isfile(path) and os.access(path, os.X_OK) else None
    for name in BINARIES:
        found = shutil.which(name)
        if found:
            return found
    return None


class AstapSolver(Solver):
    """Solve with the ASTAP command-line binary."""

    name = "astap"

    #: ASTAP solves against its own star tiles, fetched by ``astap_db``. It
    #: never looks at Astrometry.net index files, so nothing should be selected
    #: or downloaded on its behalf.
    uses_index_files = False

    @classmethod
    def available(cls, config=None):
        return find_binary(
            getattr(getattr(config, "phase2", None), "astap_path", None)
        ) is not None

    # -- solving ---------------------------------------------------------------
    def solve(self, filepath, hints, index_paths=()):
        binary = find_binary(self.config.phase2.astap_path)
        if binary is None:
            return SolveResult(
                False, backend=self.name,
                message="the ASTAP binary was not found "
                        "(apt install astap-cli, or set phase2.astap_path)",
            )

        database = self._database(hints)
        if database is None:
            return SolveResult(
                False, backend=self.name,
                message="no ASTAP star database is available for this field",
            )

        base = os.path.splitext(filepath)[0]
        self._clean(base)

        command = [binary, "-f", filepath, "-d", database, "-o", base, "-wcs"]
        command += self._hint_args(hints)

        try:
            completed = subprocess.run(
                command, capture_output=True, text=True,
                timeout=max(float(hints.timeout_s or 120.0), 10.0),
            )
        except subprocess.TimeoutExpired:
            return SolveResult(False, backend=self.name,
                               message=f"timed out after {hints.timeout_s:.0f}s")
        except OSError as exc:
            return SolveResult(False, backend=self.name, message=str(exc))

        keys = self._read_ini(base)
        if not keys.get("PLTSOLVD", "F").upper().startswith("T"):
            return SolveResult(False, backend=self.name,
                               message=self._failure_message(keys, completed))

        header = self._header(keys, filepath)
        if header is None:
            return SolveResult(False, backend=self.name,
                               message="solved, but the WCS could not be read back")

        # ASTAP tells us when the header's plate scale was wrong and what the
        # truth is. That is the same data-quality problem `solve_scale_warn_frac`
        # exists to catch, so pass it on rather than discarding it.
        warning = keys.get("WARNING")
        if warning:
            self.logger.warning("    [!] ASTAP: %s", warning)

        return SolveResult(True, header=header, matched=None,
                           backend=self.name, message=warning or None)

    # -- helpers ---------------------------------------------------------------
    def _hint_args(self, hints):
        """Translate SolveHints into ASTAP's argument conventions.

        Two conversions are easy to get wrong and silently produce a failed
        solve: ASTAP takes right ascension in **hours**, not degrees, and
        declination as **south-polar distance** (dec + 90).
        """
        args = []
        if hints.ra_deg is not None and hints.dec_deg is not None:
            args += ["-ra", f"{float(hints.ra_deg) / 15.0:.6f}",
                     "-spd", f"{float(hints.dec_deg) + 90.0:.6f}"]
        args += ["-r", f"{float(hints.radius_deg or 10.0):.3f}"]

        # Field HEIGHT in degrees. `-fov 0` asks ASTAP to work it out, which on
        # the CASSA test frames failed every time while an explicit value solved
        # in 0.1 s -- so a known scale is always passed.
        fov = 0.0
        if hints.pixel_scale and hints.naxis2:
            fov = float(hints.naxis2) * float(hints.pixel_scale) / 3600.0
        args += ["-fov", f"{fov:.5f}"]
        return args

    def _database(self, hints):
        """The tile directory to solve against, fetching what this field needs."""
        from cassa_photometry.astap_db import (
            AstapDatabaseError,
            AstapTileStore,
            resolve_db_dir,
            series_for_fov,
        )

        phase2 = self.config.phase2
        directory = resolve_db_dir(self.config)

        # Which series covers this field. "auto" derives it from the frame's
        # own height, because a database too sparse for the field is exactly
        # what produces "no solution found" on a narrow one.
        fov = None
        if hints.pixel_scale and hints.naxis2:
            fov = float(hints.naxis2) * float(hints.pixel_scale) / 3600.0
        series = (phase2.astap_db_series or "auto").strip().lower()
        if series == "auto":
            series = series_for_fov(fov) or "d50"
            self.logger.debug("    -> ASTAP series %r for a %.3f deg field.",
                              series, fov or 0.0)

        store = AstapTileStore(directory, base_url=phase2.astap_db_url,
                               download=phase2.astap_db_download,
                               logger=self.logger)
        try:
            return store.ensure(hints.ra_deg, hints.dec_deg, series=series,
                                neighbours=phase2.astap_tile_neighbours)
        except AstapDatabaseError as exc:
            self.logger.warning("    [!] %s", exc)
            # A directory the user populated themselves still works, even when
            # the incremental fetch could not run.
            return directory if os.path.isdir(directory) else None
        except Exception as exc:  # network, permissions, a moved archive
            self.logger.warning(
                "    [!] Could not top up the ASTAP database (%s). Using whatever "
                "is already in %s.", exc, directory)
            return directory if os.path.isdir(directory) else None

    def _read_ini(self, base):
        path = base + ".ini"
        keys = {}
        if not os.path.exists(path):
            return keys
        try:
            with open(path) as handle:
                for line in handle:
                    key, sep, value = line.partition("=")
                    if sep:
                        keys[key.strip()] = value.strip()
        except OSError:
            pass
        return keys

    def _failure_message(self, keys, completed):
        """Say which kind of failure this was, because the fixes differ."""
        error = keys.get("ERROR", "")
        if "star database" in error.lower():
            return ("ASTAP has no star database covering this field. The tiles "
                    "are fetched on demand; check the network, or point "
                    "phase2.astap_db_dir at a database you installed.")
        if error:
            return f"ASTAP: {error}"
        tail = (completed.stdout or completed.stderr or "").strip().splitlines()
        return tail[-1] if tail else "no solution found"

    def _header(self, keys, filepath):
        """A FITS WCS header from ASTAP's .ini solution."""
        from astropy.io import fits

        from cassa_photometry.fits_utils import open_fits

        wanted = ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
                  "CROTA1", "CROTA2", "CD1_1", "CD1_2", "CD2_1", "CD2_2")
        try:
            with open_fits(filepath) as hdul:
                header = hdul[0].header.copy()
        except Exception:
            header = fits.Header()

        found = False
        for key in wanted:
            if key in keys:
                try:
                    header[key] = float(keys[key])
                    found = True
                except ValueError:
                    continue
        if not found:
            return None

        header["CTYPE1"] = ("RA---TAN", "first parameter RA, TAN projection")
        header["CTYPE2"] = ("DEC--TAN", "second parameter DEC, TAN projection")
        header["CUNIT1"] = "deg"
        header["CUNIT2"] = "deg"
        header.setdefault("EQUINOX", 2000.0)
        return header

    def _clean(self, base):
        for suffix in TEMP_SUFFIXES:
            try:
                os.remove(base + suffix)
            except OSError:
                pass
