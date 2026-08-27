"""Data loading and robust frame statistics.

Frame statistics are computed on a central region of interest (to avoid
vignetting and amp-glow in the corners) with sigma-clipping (to reject hot /
telegraph pixels) -- the key accuracy improvement over the original full-frame
``np.mean`` / ``np.var``.
"""

import os
import glob
from collections import defaultdict

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats


def central_roi(frame, fraction=0.5):
    """Return the central ``fraction`` (per side) of a 2D array."""
    if fraction >= 1.0:
        return frame
    ny, nx = frame.shape
    hy, hx = int(ny * fraction / 2), int(nx * fraction / 2)
    cy, cx = ny // 2, nx // 2
    return frame[cy - hy:cy + hy, cx - hx:cx + hx]


def robust_mean(frame, roi_fraction=0.5, sigma=5.0):
    """Sigma-clipped mean of the central ROI (robust signal level in ADU)."""
    roi = central_roi(np.asarray(frame, dtype=np.float64), roi_fraction)
    mean, _, _ = sigma_clipped_stats(roi, sigma=sigma)
    return float(mean)


def robust_diff_variance(frame_a, frame_b, roi_fraction=0.5, sigma=5.0):
    """Signal variance from a frame pair: ``Var(A-B)/2`` on the ROI, sigma-clipped.

    Differencing removes fixed-pattern noise; sigma-clipping removes hot pixels.
    """
    diff = central_roi(np.asarray(frame_a, np.float64) - np.asarray(frame_b, np.float64), roi_fraction)
    _, _, std = sigma_clipped_stats(diff, sigma=sigma)
    return float(std ** 2) / 2.0


def load_grouped(data_dir):
    """Load calibration frames grouped by the FITS headers.

    Returns ``(bias_by_gain, flats_by_gain_exp, darks_by_temp_exp)`` where the
    values are lists of 2D arrays.
    """
    bias_by_gain = defaultdict(list)
    flats_by_gain_exp = defaultdict(lambda: defaultdict(list))
    darks_by_temp_exp = defaultdict(lambda: defaultdict(list))

    for f in glob.glob(os.path.join(data_dir, "bias", "*.fit*")):
        with fits.open(f) as hdul:
            bias_by_gain[hdul[0].header.get("GAIN", 0)].append(hdul[0].data)

    for f in glob.glob(os.path.join(data_dir, "flat", "*.fit*")):
        with fits.open(f) as hdul:
            h = hdul[0].header
            flats_by_gain_exp[h.get("GAIN", 0)][h.get("EXPTIME", 0.0)].append(hdul[0].data)

    for f in glob.glob(os.path.join(data_dir, "dark", "*.fit*")):
        with fits.open(f) as hdul:
            h = hdul[0].header
            darks_by_temp_exp[round(h.get("CCD-TEMP", 0.0))][h.get("EXPTIME", 1.0)].append(hdul[0].data)

    return bias_by_gain, flats_by_gain_exp, darks_by_temp_exp


def load_spectra(data_dir):
    """Load the QE curve and filter transmission CSVs if present.

    Returns ``(qe, filters)`` where ``qe`` is a dict or None, and ``filters`` is
    a name->dict mapping (possibly empty).
    """
    spectra_dir = os.path.join(data_dir, "spectra")
    qe = None
    qe_file = os.path.join(spectra_dir, "qe_calibration.csv")
    if os.path.exists(qe_file):
        wv, q = np.loadtxt(qe_file, delimiter=",", skiprows=1, unpack=True)
        qe = {"wavelengths": wv, "qe_percentage": q, "synthetic": False}

    filters = {}
    for f in sorted(glob.glob(os.path.join(spectra_dir, "filter_*.csv"))):
        name = os.path.basename(f).replace("filter_", "").replace(".csv", "")
        wv, unfilt, filt = np.loadtxt(f, delimiter=",", skiprows=1, unpack=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            transmission = np.where(unfilt > 0, filt / unfilt * 100.0, 0.0)
        filters[name] = {"wavelengths": wv, "transmission": transmission, "synthetic": False}

    return qe, filters
