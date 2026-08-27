"""FITS I/O for the integration phase: metadata, plane loading, master saving.

Reads the SCI/ERR/DQ planes written by phase 1 and writes the stacked master
back in the same multi-extension form.
"""

import numpy as np
from astropy.io import fits

from cassa_photometry.fits_utils import write_mef, read_mef
from cassa_photometry.instruments import ITelescopeNetworkProfile

# One shared instrument profile so filter naming matches the rest of the pipeline.
_INSTRUMENT = ITelescopeNetworkProfile()


class FITSHandler:
    """Handles FITS reading, metadata extraction, and master saving."""

    @staticmethod
    def extract_metadata(filepath):
        """Return grouping/plate-scale metadata, or None if the file is unreadable."""
        try:
            with fits.open(filepath) as hdul:
                header = hdul[0].header
                meta = {}
                meta["object"] = header.get("OBJECT", header.get("TARGET", "Unknown")).replace(" ", "").upper()
                # Canonical filter label (shared with phases 1 & 3 via the instrument profile).
                meta["filter"] = _INSTRUMENT.standardize_filter(header)
                meta["hardware_id"] = header.get("INSTRUME", "Unknown_Cam").replace(" ", "")
                raw_exp = float(header.get("EXPTIME", header.get("EXPOSURE", 0.0)))
                meta["exposure"] = int(5 * round(raw_exp / 5))
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
        bg_val = float(np.nanmedian(sci))
        header["BKGND"] = (bg_val if np.isfinite(bg_val) else 0.0, "Median background")

        write_mef(
            output_filename,
            sci=np.nan_to_num(sci, nan=0.0),
            err=None if err is None else np.nan_to_num(err, nan=0.0),
            dq=dq,
            header=header,
            history="PHASE 2: SNR anchor, rejection, inverse-variance weighted; SCI/ERR/DQ",
        )
        return pixel_scale
