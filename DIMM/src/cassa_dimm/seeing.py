"""DIMM physics: differential image motion -> Fried parameter -> seeing.

Implements the standard Sarazin & Roddier (1990) / Tokovinin (2002) reduction,
generalised so the longitudinal axis follows an arbitrary prism-doublet vector,
with airmass correction to the zenith and an SNR-weighted combination of many
star doublets.
"""

from dataclasses import dataclass, field

import numpy as np
from astropy.stats import sigma_clip

_ARCSEC_TO_RAD = 4.84813681109536e-6
_SEEING_CONST = 0.98            # FWHM = 0.98 * lambda / r0


# --- Response coefficients --------------------------------------------------- #
def response_coefficients(aperture_diameter_m, hole_separation_m):
    """Sarazin & Roddier (1990) differential-variance coefficients ``(c_l, c_t)``.

    The differential image-motion variance is
    ``sigma^2 = 2 * lambda^2 * r0^(-5/3) * c`` with

    * ``c_l = 0.179 D^(-1/3) - 0.0968 d^(-1/3)`` (longitudinal, along the baseline)
    * ``c_t = 0.179 D^(-1/3) - 0.1450 d^(-1/3)`` (transverse)
    """
    D = aperture_diameter_m
    d = hole_separation_m
    c_l = 0.179 * D ** (-1 / 3) - 0.0968 * d ** (-1 / 3)
    c_t = 0.179 * D ** (-1 / 3) - 0.1450 * d ** (-1 / 3)
    return c_l, c_t


def variance_to_r0(variance_rad2, coefficient, wavelength_m):
    """Convert a differential-motion variance (rad^2) to the Fried parameter r0 (m)."""
    if variance_rad2 <= 0 or coefficient <= 0:
        return np.nan
    return (2.0 * wavelength_m ** 2 * coefficient / variance_rad2) ** 0.6


def r0_to_seeing_arcsec(r0_m, wavelength_m):
    """Seeing FWHM in arcseconds from r0 (at the given wavelength)."""
    if not np.isfinite(r0_m) or r0_m <= 0:
        return np.nan
    return _SEEING_CONST * (wavelength_m / r0_m) / _ARCSEC_TO_RAD


# --- Airmass ----------------------------------------------------------------- #
def kasten_young_airmass(altitude_deg):
    """Airmass via Kasten & Young (1989); robust down to the horizon."""
    if altitude_deg is None or altitude_deg <= 0:
        return np.nan
    z = 90.0 - altitude_deg
    return 1.0 / (np.cos(np.radians(z)) + 0.50572 * (altitude_deg + 6.07995) ** -1.6364)


def correct_seeing_to_zenith(seeing_arcsec, airmass):
    """Scale an observed seeing to the zenith: ``eps0 = eps * X^(-0.6)``."""
    if not np.isfinite(seeing_arcsec) or airmass is None or not np.isfinite(airmass) or airmass <= 0:
        return np.nan
    return seeing_arcsec * airmass ** (-0.6)


# --- Per-doublet reduction --------------------------------------------------- #
@dataclass
class DoubletResult:
    """Seeing measured from one star's two prism spots."""

    seeing_l: float
    seeing_t: float
    seeing: float           # mean of longitudinal & transverse
    r0_l: float
    r0_t: float
    var_l_rad2: float
    var_t_rad2: float
    n_used: int
    snr: float
    seeing_err: float
    flags: list = field(default_factory=list)


def measure_doublet(dx_pix, dy_pix, prism_unit_vec, hardware, qc, snr=1.0):
    """Reduce one doublet's per-frame separation offsets to a seeing estimate.

    Parameters
    ----------
    dx_pix, dy_pix : ndarray
        Per-frame separation ``spot1 - spot2`` in pixels.
    prism_unit_vec : (float, float)
        Unit vector along the doublet (defines the longitudinal axis).
    hardware : HardwareConfig
    qc : QCConfig
    snr : float
        Mean SNR of the two spots (used as the combination weight).
    """
    dx = np.asarray(dx_pix, dtype=float)
    dy = np.asarray(dy_pix, dtype=float)
    ux, uy = prism_unit_vec

    # Project onto longitudinal (along baseline) and transverse axes.
    long_comp = dx * ux + dy * uy
    tran_comp = -dx * uy + dy * ux

    # Reject outlier frames (clouds, glitches) independently on each axis.
    long_c = sigma_clip(long_comp, sigma=qc.sigma_clip, maxiters=3, masked=True)
    tran_c = sigma_clip(tran_comp, sigma=qc.sigma_clip, maxiters=3, masked=True)
    n_used = int(min(long_c.count(), tran_c.count()))

    scale2 = hardware.scale_rad ** 2
    var_l = float(np.ma.var(long_c)) * scale2
    var_t = float(np.ma.var(tran_c)) * scale2

    c_l, c_t = response_coefficients(hardware.aperture_diameter_m, hardware.hole_separation_m)
    r0_l = variance_to_r0(var_l, c_l, hardware.wavelength_m)
    r0_t = variance_to_r0(var_t, c_t, hardware.wavelength_m)
    seeing_l = r0_to_seeing_arcsec(r0_l, hardware.wavelength_m)
    seeing_t = r0_to_seeing_arcsec(r0_t, hardware.wavelength_m)

    seeing = float(np.nanmean([seeing_l, seeing_t]))
    # Relative seeing error: seeing ~ (sigma^2)^{3/5}, so d(eps)/eps = 0.6*sqrt(2/N).
    seeing_err = seeing * 0.6 * np.sqrt(2.0 / max(n_used, 1)) if np.isfinite(seeing) else np.nan

    flags = []
    if np.isfinite(seeing_l) and np.isfinite(seeing_t) and seeing > 0:
        if abs(seeing_l - seeing_t) / seeing > qc.max_lt_discrepancy:
            flags.append("LT_DISCREPANCY")
    return DoubletResult(seeing_l, seeing_t, seeing, r0_l, r0_t, var_l, var_t,
                         n_used, float(snr), seeing_err, flags)


# --- Window combination ------------------------------------------------------ #
@dataclass
class WindowResult:
    """Combined seeing from all valid doublets in a window."""

    seeing_zenith: float
    seeing_raw: float
    seeing_l: float
    seeing_t: float
    r0_cm: float
    airmass: float
    n_frames: int
    n_stars: int
    star_scatter: float
    mean_snr: float
    seeing_err: float
    flags: list = field(default_factory=list)


def combine_window(doublets, airmass, n_frames, qc):
    """SNR-weighted combination of per-doublet results into a window result."""
    good = [d for d in doublets if np.isfinite(d.seeing) and d.n_used >= qc.min_frames]
    flags = []
    if len(good) < qc.min_doublets:
        flags.append("TOO_FEW_DOUBLETS")
        return WindowResult(np.nan, np.nan, np.nan, np.nan, np.nan,
                            airmass if airmass else np.nan, n_frames, len(good),
                            np.nan, np.nan, np.nan, flags)

    w = np.array([max(d.snr, 1e-3) for d in good])
    seeings = np.array([d.seeing for d in good])
    seeing_l = float(np.average([d.seeing_l for d in good], weights=w))
    seeing_t = float(np.average([d.seeing_t for d in good], weights=w))
    seeing_raw = float(np.average(seeings, weights=w))
    r0_cm = float(np.average([0.5 * (d.r0_l + d.r0_t) for d in good], weights=w)) * 100.0
    mean_snr = float(np.mean(w))

    # Uncertainty: star-to-star scatter combined with the mean per-doublet error.
    star_scatter = float(np.std(seeings, ddof=1)) if len(good) > 1 else np.nan
    per_doublet_err = float(np.mean([d.seeing_err for d in good]))
    if len(good) > 1:
        seeing_err = float(np.hypot(star_scatter / np.sqrt(len(good)), per_doublet_err))
    else:
        seeing_err = per_doublet_err

    seeing_zenith = correct_seeing_to_zenith(seeing_raw, airmass)
    if not np.isfinite(seeing_zenith):
        seeing_zenith = seeing_raw
        flags.append("NO_AIRMASS")
    if any("LT_DISCREPANCY" in d.flags for d in good):
        flags.append("LT_DISCREPANCY")

    return WindowResult(seeing_zenith, seeing_raw, seeing_l, seeing_t, r0_cm,
                        airmass if airmass else np.nan, n_frames, len(good),
                        star_scatter, mean_snr, seeing_err, flags)
