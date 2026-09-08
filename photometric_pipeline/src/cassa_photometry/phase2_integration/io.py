"""FITS I/O for the integration phase: metadata, plane loading, master saving.

Reads the SCI/ERR/DQ planes written by phase 1 and writes the stacked master
back in the same multi-extension form.
"""

import numpy as np
from astropy.io import fits

from cassa_photometry.fits_utils import CALVERS, read_mef, write_mef
from cassa_photometry.instruments import get_profile


class FITSHandler:
    """Handles FITS reading, metadata extraction, and master saving."""

    @staticmethod
    def extract_metadata(filepath, instrument=None, config=None):
        """Return grouping/plate-scale metadata, or None if the file is unreadable.

        ``instrument`` supplies the canonical filter label, shared with phases 1
        and 3 so the three phases group frames identically. Pass ``config``
        instead to resolve the profile from a configuration; passing **neither**
        gives the ``generic`` profile, which is a different thing from "the
        profile this run was configured with" -- every caller inside the
        pipeline passes one.
        """
        instrument = instrument or get_profile(
            getattr(config, "instrument", None), config=config
        )
        try:
            with fits.open(filepath) as hdul:
                header = hdul[0].header
                meta = {}
                meta["object"] = header.get("OBJECT", header.get("TARGET", "Unknown")).replace(" ", "").upper()
                meta["filter"] = instrument.standardize_filter(header)
                meta["hardware_id"] = header.get("INSTRUME", "Unknown_Cam").replace(" ", "")
                raw_exp = float(header.get("EXPTIME", header.get("EXPOSURE", 0.0)))
                # Quantised only for *grouping*, so 118 s and 120 s frames stack
                # together. The true value is kept: normalising by the quantised
                # one would put a ~2% (0.02 mag) systematic straight into every
                # magnitude.
                meta["exposure"] = int(5 * round(raw_exp / 5))
                meta["exposure_true"] = raw_exp
                meta["date_obs"] = str(header.get("DATE-OBS", "") or "")
                meta["mjd_obs"] = header.get("MJD-OBS")
                meta["site_long"] = header.get("SITELONG")
                meta["airmass"] = header.get("AIRMASS")
                meta["pixel_size_um"] = float(header.get("XPIXSZ", 0.0))
                meta["focal_length_mm"] = float(header.get("FOCALLEN", 0.0))
                meta["binning"] = int(header.get("XBINNING", 1))
                return meta
        except Exception:
            return None

    @staticmethod
    def load_data(filepath):
        """Return only the SCI plane as float32 (for detection / anchor selection)."""
        sci, _, _, _ = read_mef(filepath)
        return np.asarray(sci, dtype=np.float32)

    @staticmethod
    def load_planes(filepath):
        """Return ``(sci, err, dq)`` as float32/int arrays; ``err``/``dq`` may be None."""
        sci, err, dq, _ = read_mef(filepath)
        sci = np.asarray(sci, dtype=np.float32)
        err = None if err is None else np.asarray(err, dtype=np.float32)
        dq = None if dq is None else np.asarray(dq, dtype=np.int32)
        return sci, err, dq

    @staticmethod
    def save_master(sci, err, dq, ref_filepath, output_filename, num_frames, meta):
        """Write the stacked master as a SCI/ERR/DQ multi-extension FITS file."""
        with fits.open(ref_filepath) as hdul:
            header = hdul[0].header.copy()

        pixel_scale = 0.0
        if meta["focal_length_mm"] > 0 and meta["pixel_size_um"] > 0:
            pixel_scale = (meta["pixel_size_um"] / meta["focal_length_mm"]) * 206.265

        header["STACKCNT"] = (num_frames, "Total frames integrated")
        header["TOT_EXP"] = (num_frames * meta["exposure"], "[s] Total integration time")
        header["PIXSCALE"] = (round(pixel_scale, 4), "[arcsec/pix] Plate scale")
        header["BUNIT"] = "electron"
        header["CALVERS"] = (CALVERS, "Calibration vintage")
        bg_val = float(np.nanmedian(sci))
        header["BKGND"] = (bg_val if np.isfinite(bg_val) else 0.0, "Median background")

        # The stack is a weighted MEAN, so a pixel holds electrons per single
        # frame -- not per TOT_EXP. Normalising by TOT_EXP would be wrong by a
        # factor of N. EXPMEAN is the true mean single-frame exposure, unrounded,
        # and is the number to divide by.
        stats = meta.get("stack_stats") or {}
        if stats.get("exposure_mean"):
            header["EXPMEAN"] = (round(float(stats["exposure_mean"]), 4),
                                 "[s] Mean single-frame exposure (divide by THIS)")
        if stats.get("airmass_mean"):
            header["AIRMASS"] = (round(float(stats["airmass_mean"]), 4),
                                 "Mean airmass of the contributing frames")
        if stats.get("background_level") is not None:
            header["BKGLEVEL"] = (round(float(stats["background_level"]), 4),
                                  "[e-] Sky level removed before stacking")
        for card, key, comment in (
            ("FWHMPX", "fwhm", "[pixel] Measured PSF FWHM of this master"),
            ("FWHMSTD", "fwhm_std", "[pixel] Frame-to-frame FWHM scatter"),
            ("FWHMMIN", "fwhm_min", "[pixel] Best contributing frame FWHM"),
            ("FWHMMAX", "fwhm_max", "[pixel] Worst contributing frame FWHM"),
            ("NFWHMREJ", "n_fwhm_rejected", "Frames rejected for poor seeing"),
        ):
            if stats.get(key) is not None:
                header[card] = (stats[key], comment)
            else:
                # The header was copied from the anchor *frame*, which phase 1
                # stamped with its own FWHMPX. Leaving that in place would hand
                # phase 3 one frame's seeing as if it were the stack's, and
                # phase 3 sizes every aperture from this card -- so an
                # unmeasured PSF must read as absent, not as the anchor's.
                header.remove(card, ignore_missing=True)
        if stats.get("epoch"):
            header["EPOCHKEY"] = (stats["epoch"], "Observing-night bin this master covers")

        # Say what was not done, in the file itself -- the same contract phase 1
        # honours with CALSKIP. A stack that skipped background subtraction or
        # registration must not be mistakable for a default one later.
        if meta.get("steps_skipped"):
            header["STEPSKIP"] = (",".join(meta["steps_skipped"]),
                                  "Phase 2 steps deliberately skipped")
        if meta.get("steps_plan"):
            header["STEPPLAN"] = (",".join(meta["steps_plan"]),
                                  "Phase 2 steps run, in order")

        write_mef(
            output_filename,
            # SCI is zero-filled where there is no coverage so that display and
            # source detection behave; the DQ NO_DATA bit is what says those
            # pixels are not measurements.
            sci=np.nan_to_num(sci, nan=0.0),
            # ERR is NOT zero-filled. A zero uncertainty is infinite
            # inverse-variance weight to everything downstream, so an uncovered
            # pixel would be trusted more than any real measurement. NaN is the
            # honest value, and every consumer already has to handle it.
            err=err,
            dq=dq,
            header=header,
            history="PHASE 2: SNR anchor, rejection, inverse-variance weighted; SCI/ERR/DQ",
        )
        return pixel_scale
