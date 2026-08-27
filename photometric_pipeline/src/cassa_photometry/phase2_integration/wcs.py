"""Astrometric calibration via a local Astrometry.net ``solve-field`` install.

``solve-field`` only understands a single-image FITS, and its ``--new-fits``
output drops extra extensions. So we solve on the SCI (primary) plane, then
re-attach the ERR/DQ planes to the WCS-solved result, keeping the master a full
SCI/ERR/DQ multi-extension file.
"""

import os
import shutil
import tempfile
import subprocess

from astropy.io import fits

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import read_mef, write_mef


class WCSSolver:
    def __init__(self, logger, config=None):
        self.logger = logger
        self.config = config or load_config()
        self.index_dir = self.config.resolve_astrometry_index_dir()
        self.global_anchor_ra = None
        self.global_anchor_dec = None
        self.config_path = self._write_config()

        if not os.path.isdir(self.index_dir):
            self.logger.warning(
                f"[!] Astrometry index directory not found: {self.index_dir}. "
                "Set CASSA_ASTROMETRY_INDEX or phase2.astrometry_index_dir. WCS solving will fail."
            )

    def _write_config(self):
        """Write a fresh solve-field config pointing at the resolved index dir."""
        fd, path = tempfile.mkstemp(prefix="astrometry_", suffix=".cfg")
        with os.fdopen(fd, "w") as handle:
            handle.write(f"add_path {self.index_dir}\nautoindex\ninparallel\n")
        return path

    def solve(self, group, is_rescue=False, keep_temps=False):
        filepath = group.master_filepath
        if not filepath or not os.path.exists(filepath):
            return False

        # Preserve the error/quality planes across the solve-field round-trip.
        _, err, dq, _ = read_mef(filepath)

        output_fits = os.path.splitext(filepath)[0] + "_wcs.fits"
        flag_file = os.path.splitext(filepath)[0] + ".solved"
        for f in (output_fits, flag_file):
            if os.path.exists(f):
                os.remove(f)

        scale_args, pix_scale = [], 0.0
        try:
            with fits.open(filepath) as hdul:
                hdr = hdul[0].header
                pix_scale = float(hdr.get("PIXSCALE", 0.0))
                ra_hint = self.global_anchor_ra or hdr.get("OBJCTRA") or hdr.get("RA")
                dec_hint = self.global_anchor_dec or hdr.get("OBJCTDEC") or hdr.get("DEC")
            if pix_scale > 0:
                scale_args += ["--scale-units", "arcsecperpix",
                               "--scale-low", str(pix_scale * 0.75), "--scale-high", str(pix_scale * 1.25)]
            if ra_hint and dec_hint:
                ra_f = str(ra_hint).strip().replace(" ", ":")
                dec_f = str(dec_hint).strip().replace(" ", ":")
                scale_args += ["--ra", ra_f, "--dec", dec_f, "--radius", "3.0"]
        except Exception:
            pass

        if pix_scale > 1.0:
            ds, obj, sig = "4", "150", "10"
        elif 0.0 < pix_scale <= 1.0:
            if any(nb in filepath for nb in ["_Ha_", "_OIII_", "_SII_"]):
                ds, obj, sig = "1", "1000", "3"
            else:
                ds, obj, sig = "2", "600", "8"
        else:
            ds, obj, sig = "2", "300", "8"

        self.logger.info(f"    -> WCS Profile: DS {ds}x | Objs: {obj} | Sigma: {sig}")
        base_cmd = ["solve-field", filepath, "--config", self.config_path, "--overwrite",
                    "--no-plots", "--downsample", ds, "--objs", obj, "--sigma", sig,
                    "--new-fits", output_fits]

        try:
            res_fast = subprocess.run(base_cmd + ["--cpulimit", "60"] + scale_args,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            solved = os.path.exists(flag_file) and os.path.exists(output_fits)

            if not solved and not is_rescue:
                self.logger.info("    [!] Fast solve failed. Engaging Blind Solve (300s)...")
                self._log_solve_output("FAST", res_fast)
                res_blind = subprocess.run(base_cmd + ["--cpulimit", "300"],
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                solved = os.path.exists(flag_file) and os.path.exists(output_fits)
                if not solved:
                    self._log_solve_output("BLIND", res_blind)
        except FileNotFoundError:
            self.logger.error("    [-] FATAL WCS ERROR: 'solve-field' not found. Install Astrometry.net.")
            return False

        if os.path.exists(flag_file) and os.path.exists(output_fits):
            # Re-attach ERR/DQ to the WCS-solved SCI plane.
            with fits.open(output_fits) as hdul:
                solved_data = hdul[0].data
                solved_header = hdul[0].header.copy()
            self.global_anchor_ra = solved_header.get("CRVAL1")
            self.global_anchor_dec = solved_header.get("CRVAL2")
            write_mef(filepath, sci=solved_data, err=err, dq=dq, header=solved_header,
                      history="PHASE 2: WCS solved via Astrometry.net")
            os.remove(output_fits)
            if not keep_temps:
                self._clean_temps(filepath)
            group.wcs_solved = True
            self.logger.info("    [+] SUCCESS: WCS mapped.")
            return True

        self.logger.error("    [-] WCS FAILED.")
        if not keep_temps:
            self._clean_temps(filepath)
        return False

    def _log_solve_output(self, label, result):
        self.logger.error(f"        --- ASTROMETRY {label} SOLVE LOG ---")
        if result.stdout:
            self.logger.error(result.stdout.strip())
        if result.stderr:
            self.logger.error(result.stderr.strip())
        self.logger.error("        ---------------------------------")

    def _clean_temps(self, base_filepath):
        base = os.path.splitext(base_filepath)[0]
        for ext in ["-indx.xyls", ".axy", ".corr", ".match", ".rdls", ".solved", ".wcs"]:
            try:
                os.remove(f"{base}{ext}")
            except OSError:
                pass
