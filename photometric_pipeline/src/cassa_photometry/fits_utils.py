"""Multi-extension FITS (MEF) I/O and the data-quality (DQ) bit conventions.

Every science image the pipeline produces carries three planes, following the
convention used by mature reduction pipelines (Astropy ``ccdproc``, Gemini
DRAGONS, LCO BANZAI):

* ``SCI`` -- the science data, held in the **primary** HDU (so ordinary tools
  such as ``fits.getdata`` and ``solve-field`` keep working).
* ``ERR`` -- the 1-sigma uncertainty, in the same units as ``SCI``.
* ``DQ``  -- an integer bitmask flagging bad pixels (see the bit constants).
"""

import numpy as np
from astropy.io import fits

# --- Data-quality bit flags ---------------------------------------------------
DQ_GOOD = 0
DQ_SATURATED = 1      # pixel at/above the saturation level
DQ_BAD_PIXEL = 2      # hot/dead pixel from the bad-pixel mask
DQ_COSMIC_RAY = 4     # flagged by cosmic-ray rejection
DQ_NO_DATA = 8        # NaN / no coverage (e.g. outside the warp footprint)
DQ_REJECTED = 16      # some contributing frames were rejected here (stacks only)

DQ_FLAG_NAMES = {
    DQ_SATURATED: "SATURATED",
    DQ_BAD_PIXEL: "BAD_PIXEL",
    DQ_COSMIC_RAY: "COSMIC_RAY",
    DQ_NO_DATA: "NO_DATA",
    DQ_REJECTED: "REJECTED",
}

ERR_EXTNAME = "ERR"
DQ_EXTNAME = "DQ"

# --- Calibration vintage ------------------------------------------------------
# Stamped as CALVERS on every product the pipeline writes. Bump it whenever a
# change makes new output numerically incomparable with old output, so a file
# states its own calibration vintage and mixed vintages are detectable rather
# than silently averaged together. Vintage 1 is anything written before this
# card existed (no CALVERS => 1).
#
#   1  original release
#   2  scientific hardening: DQ propagated through stacking, aperture-corrected
#      total-flux zero point, per-band photometric system, exposure
#      normalisation. See docs/CHANGELOG.md and docs/master-plan.md.
CALVERS = 2


def build_dq(shape, saturated=None, bad_pixel=None, cosmic_ray=None, no_data=None,
             rejected=None):
    """Combine boolean masks into a single integer DQ bitmask array."""
    dq = np.zeros(shape, dtype=np.int32)
    for mask, flag in (
        (saturated, DQ_SATURATED),
        (bad_pixel, DQ_BAD_PIXEL),
        (cosmic_ray, DQ_COSMIC_RAY),
        (no_data, DQ_NO_DATA),
        (rejected, DQ_REJECTED),
    ):
        if mask is not None:
            dq[np.asarray(mask, dtype=bool)] |= flag
    return dq


def write_mef(path, sci, err=None, dq=None, header=None, overwrite=True, history=None):
    """Write a SCI/ERR/DQ multi-extension FITS file.

    Parameters
    ----------
    path : str
        Output filename.
    sci : ndarray
        Science plane (stored in the primary HDU).
    err : ndarray, optional
        1-sigma uncertainty plane.
    dq : ndarray, optional
        Integer data-quality bitmask.
    header : astropy.io.fits.Header, optional
        Header for the primary HDU. WCS and instrument keywords should live here.
    history : str or list of str, optional
        HISTORY card(s) to append.
    """
    hdr = header.copy() if header is not None else fits.Header()
    if history:
        for line in ([history] if isinstance(history, str) else history):
            hdr["HISTORY"] = line

    hdus = [fits.PrimaryHDU(data=np.asarray(sci, dtype=np.float32), header=hdr)]
    if err is not None:
        hdus.append(fits.ImageHDU(data=np.asarray(err, dtype=np.float32), name=ERR_EXTNAME))
    if dq is not None:
        hdus.append(fits.ImageHDU(data=np.asarray(dq, dtype=np.int32), name=DQ_EXTNAME))

    fits.HDUList(hdus).writeto(path, overwrite=overwrite)


def read_mef(path):
    """Read a SCI/ERR/DQ FITS file.

    Returns
    -------
    sci : ndarray
    err : ndarray or None
    dq : ndarray or None
    header : astropy.io.fits.Header
        The primary (SCI) header.

    Notes
    -----
    Tolerant of legacy single-HDU files: ``err`` and ``dq`` come back ``None``
    when the extensions are absent so downstream code can fall back gracefully.
    """
    with fits.open(path) as hdul:
        sci = hdul[0].data
        header = hdul[0].header
        err = hdul[ERR_EXTNAME].data if ERR_EXTNAME in hdul else None
        dq = hdul[DQ_EXTNAME].data if DQ_EXTNAME in hdul else None
        # Realise the arrays before the file closes.
        sci = None if sci is None else np.asarray(sci, dtype=np.float32)
        err = None if err is None else np.asarray(err, dtype=np.float32)
        dq = None if dq is None else np.asarray(dq, dtype=np.int32)
    return sci, err, dq, header
