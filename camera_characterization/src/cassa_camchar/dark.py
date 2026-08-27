"""Dark current (temperature sweep) and dark linearity / amp-glow analysis."""

import numpy as np

from cassa_camchar.io import robust_mean


def dark_current(dark_frame, bias_mean, exposure_s, gain, config):
    """Dark current in e-/pixel/s: ``(mean(dark) - bias) * gain / t``."""
    a = config.analysis
    signal_e = (robust_mean(dark_frame, a.roi_fraction, a.sigma_clip) - bias_mean) * gain
    return signal_e / exposure_s if exposure_s > 0 else np.nan


def thermal_electrons(dark_frame, bias_mean, gain, config):
    """Total accumulated thermal signal in electrons."""
    a = config.analysis
    return (robust_mean(dark_frame, a.roi_fraction, a.sigma_clip) - bias_mean) * gain


def temperature_sweep(darks_by_temp_exp, bias_mean, gain, config):
    """Return ``(temperatures, dark_currents)`` using the longest exposure per temperature."""
    temps, currents = [], []
    for t in sorted(darks_by_temp_exp.keys()):
        longest = max(darks_by_temp_exp[t].keys())
        frame = darks_by_temp_exp[t][longest][0]
        temps.append(t)
        currents.append(dark_current(frame, bias_mean, longest, gain, config))
    return np.array(temps, dtype=float), np.array(currents, dtype=float)


def linearity(darks_by_temp_exp, bias_mean, gain, config):
    """Fit total thermal electrons vs exposure time at the linearity temperature.

    Returns a dict with the fit and a non-linearity residual (%) as an amp-glow
    indicator, or ``None`` if there are fewer than two exposures (no fabricated
    fallback, unlike the original script).
    """
    if not darks_by_temp_exp:
        return None
    target = min(darks_by_temp_exp.keys(), key=lambda k: abs(k - config.capture.linearity_temp_c))
    exps, elec = [], []
    for exp in sorted(darks_by_temp_exp[target].keys()):
        exps.append(exp)
        elec.append(thermal_electrons(darks_by_temp_exp[target][exp][0], bias_mean, gain, config))

    if len(exps) < 2:
        return None
    exps, elec = np.array(exps, dtype=float), np.array(elec, dtype=float)
    slope, intercept = np.polyfit(exps, elec, 1)
    residual = elec - (slope * exps + intercept)
    rng = float(elec.max() - elec.min())
    nonlinearity_pct = float(np.max(np.abs(residual)) / rng * 100.0) if rng > 0 else 0.0
    return {"exposures": exps, "thermal_electrons": elec, "temperature_c": float(target),
            "slope": float(slope), "intercept": float(intercept),
            "nonlinearity_pct": nonlinearity_pct}
