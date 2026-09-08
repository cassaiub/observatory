"""The Astrometry.net ``solve-field`` backend.

Preferred when the binary is installed: it is mature, and its ``.corr`` output
gives the matched star pairs the astrometric residual is computed from.

Three things here are deliberate departures from how this was previously
invoked, each of them a bug that was costing solves:

* **The config file is written per solve and deleted afterwards.** The old code
  left one ``/tmp/astrometry_*.cfg`` behind per solver instance.
* **Index files are named explicitly** rather than pointing the solver at a
  directory with ``autoindex``. That is what makes on-demand fetching possible,
  and it also stops a user with the full 5 GB set from loading all of it.
* **The noise level is measured, not assumed.** ``--sigma`` is not a threshold
  multiplier -- it is the assumed image noise in ADU, so a hardcoded value is
  detector-specific. On CASSA-like frames, assuming 5 ADU where the truth is 56
  made the extractor return 499 sources of which 28 were real, burying the stars
  it needed among noise peaks.
"""

import os
import shutil
import subprocess
import tempfile

import numpy as np
from astropy.io import fits

from cassa_photometry.phase2_integration.solvers.base import Solver, SolveResult

#: Scratch files solve-field leaves beside the frame.
TEMP_SUFFIXES = ("-indx.xyls", ".axy", ".corr", ".match", ".rdls", ".solved",
                 ".wcs", ".new")


class SolveFieldSolver(Solver):
    """Shell out to ``solve-field``."""

    name = "solve-field"

    @classmethod
    def available(cls):
        return shutil.which("solve-field") is not None

    def solve(self, filepath, hints, index_paths=()):
        base = os.path.splitext(filepath)[0]
        output_fits = base + "_wcs.fits"
        flag_file = base + ".solved"
        for path in (output_fits, flag_file):
            if os.path.exists(path):
                os.remove(path)

        config_path = self._write_config(index_paths)
        try:
            command = self._command(filepath, output_fits, hints, config_path)
            self.logger.debug("solve-field: %s", " ".join(command))
            try:
                completed = subprocess.run(
                    command, capture_output=True, text=True,
                    timeout=max(hints.timeout_s * 1.5, 30.0),
                )
            except FileNotFoundError:
                return SolveResult(False, backend=self.name,
                                   message="'solve-field' is not installed")
            except subprocess.TimeoutExpired:
                return SolveResult(False, backend=self.name, message="solve-field timed out")

            solved = os.path.exists(flag_file) and os.path.exists(output_fits)
            if not solved:
                return SolveResult(
                    False, backend=self.name,
                    message=_last_lines(completed.stdout, completed.stderr),
                )

            with fits.open(output_fits) as hdul:
                header = hdul[0].header.copy()
            os.remove(output_fits)
            return SolveResult(True, header=header, matched=self._corr(base),
                               backend=self.name)
        finally:
            try:
                os.remove(config_path)
            except OSError:
                pass

    # -- helpers ---------------------------------------------------------------
    def _write_config(self, index_paths):
        """A config naming exactly the index files this solve should consider."""
        handle, path = tempfile.mkstemp(prefix="cassa_astrometry_", suffix=".cfg")
        with os.fdopen(handle, "w") as fh:
            fh.write("inparallel\n")
            if index_paths:
                for index_path in index_paths:
                    fh.write(f"index {index_path}\n")
            else:
                # Nothing was selected: fall back to whatever the configured
                # directory holds, which is exactly the old behaviour.
                fh.write(f"add_path {self.config.resolve_astrometry_index_dir()}\n")
                fh.write("autoindex\n")
        return path

    def _command(self, filepath, output_fits, hints, config_path):
        downsample = "2" if (hints.naxis1 or 0) > 1500 else "1"
        command = [
            "solve-field", filepath, "--config", config_path, "--overwrite",
            "--no-plots", "--downsample", downsample, "--objs", "200",
            "--new-fits", output_fits, "--cpulimit", str(int(hints.timeout_s)),
        ]
        if hints.noise_adu and hints.noise_adu > 0:
            # Measured, not assumed. See the module docstring.
            command += ["--sigma", f"{hints.noise_adu:.4g}"]
        if hints.scale_low and hints.scale_high:
            command += ["--scale-units", "arcsecperpix",
                        "--scale-low", f"{hints.scale_low:.6g}",
                        "--scale-high", f"{hints.scale_high:.6g}"]
        if hints.ra_deg is not None and hints.dec_deg is not None:
            command += ["--ra", f"{float(hints.ra_deg):.6f}",
                        "--dec", f"{float(hints.dec_deg):.6f}",
                        "--radius", f"{float(hints.radius_deg):.3f}"]
        return command

    def _corr(self, base):
        """Matched field/index star pairs from the ``.corr`` table."""
        corr_path = base + ".corr"
        if not os.path.exists(corr_path):
            return None
        try:
            corr = fits.getdata(corr_path)
            return (np.asarray(corr["field_ra"], dtype=float),
                    np.asarray(corr["field_dec"], dtype=float),
                    np.asarray(corr["index_ra"], dtype=float),
                    np.asarray(corr["index_dec"], dtype=float))
        except Exception as exc:
            self.logger.warning("Could not read %s: %s", os.path.basename(corr_path), exc)
            return None

    @staticmethod
    def clean_temps(filepath):
        base = os.path.splitext(filepath)[0]
        for suffix in TEMP_SUFFIXES:
            try:
                os.remove(base + suffix)
            except OSError:
                pass


def _last_lines(stdout, stderr, n=3):
    """The tail of the solver's output, for a one-line failure message."""
    text = (stderr or "").strip() or (stdout or "").strip()
    lines = [line for line in text.splitlines() if line.strip()]
    return " | ".join(lines[-n:]) if lines else "no output"
