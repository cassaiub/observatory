"""Photon Transfer Curve: system gain, read noise, full well, dynamic range."""

import numpy as np

from cassa_camchar.io import robust_mean, robust_diff_variance


def analyze_ptc(flat_pairs, config):
    """Analyse one gain setting's flat pairs.

    Parameters
    ----------
    flat_pairs : dict
        ``exptime -> list of frames`` (each exposure needs >= 2 frames).
    config : CamCharConfig

    Returns
    -------
    dict with ``means``, ``variances`` (ADU / ADU^2), ``gain`` (e-/ADU),
    ``full_well_adu``, ``full_well_e``, and the boolean ``linear_mask``.
    """
    a = config.analysis
    means, variances = [], []
    for _exp, frames in flat_pairs.items():
        if len(frames) < 2:
            continue
        m = 0.5 * (robust_mean(frames[0], a.roi_fraction, a.sigma_clip)
                   + robust_mean(frames[1], a.roi_fraction, a.sigma_clip))
        v = robust_diff_variance(frames[0], frames[1], a.roi_fraction, a.sigma_clip)
        means.append(m)
        variances.append(v)

    means, variances = np.array(means), np.array(variances)
    if means.size < 2:
        return {"means": means, "variances": variances, "gain": np.nan,
                "full_well_adu": np.nan, "full_well_e": np.nan,
                "linear_mask": np.zeros_like(means, dtype=bool)}

    order = np.argsort(means)
    means, variances = means[order], variances[order]

    # System gain = 1 / slope of variance-vs-mean over the linear region.
    mx = float(means.max())
    linear_mask = (means > a.ptc_linear_lo * mx) & (means < a.ptc_linear_hi * mx)
    if linear_mask.sum() >= 2:
        slope, _ = np.polyfit(means[linear_mask], variances[linear_mask], 1)
        gain = 1.0 / slope if slope > 0 else np.nan
    else:
        gain = np.nan

    full_well_adu = _full_well_adu(means, variances)
    full_well_e = full_well_adu * gain if np.isfinite(gain) else np.nan
    return {"means": means, "variances": variances, "gain": gain,
            "full_well_adu": full_well_adu, "full_well_e": full_well_e,
            "linear_mask": linear_mask}


def read_noise_electrons(bias_a, bias_b, gain, config):
    """Read noise (e-) from a bias pair: ``std(A-B)/sqrt(2) * gain`` (robust)."""
    a = config.analysis
    var_signal = robust_diff_variance(bias_a, bias_b, a.roi_fraction, a.sigma_clip)
    return float(np.sqrt(max(var_signal, 0.0)) * gain)


def dynamic_range_stops(full_well_e, read_noise_e):
    """Dynamic range in stops: ``log2(FWC / RN)``."""
    if not (np.isfinite(full_well_e) and read_noise_e > 0):
        return np.nan
    return float(np.log2(full_well_e / read_noise_e))


def _full_well_adu(means, variances):
    """Signal (ADU) at the PTC variance turnover, lightly smoothed for robustness."""
    if variances.size < 3:
        return float(means[np.argmax(variances)]) if variances.size else np.nan
    smooth = np.convolve(variances, np.ones(3) / 3.0, mode="same")
    return float(means[np.argmax(smooth)])
