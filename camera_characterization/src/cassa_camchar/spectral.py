"""Quantum efficiency and filter transmission analysis.

Real lab CSVs (produced by a monochromator sweep) are used when present. When
they are absent the pipeline falls back to a documented theoretical model, but
tags the result ``synthetic=True`` so it is never mistaken for a measurement.
"""

import numpy as np

from cassa_camchar.io import load_spectra

# Theoretical filter profiles (center nm, width nm, power, peak fraction) for the
# LRGB + narrowband set. Shared with the simulator so both agree.
FILTER_PROFILES = {
    "Luminance": (550.0, 150.0, 6, 0.98),
    "Red": (640.0, 55.0, 6, 0.97),
    "Green": (535.0, 45.0, 6, 0.96),
    "Blue": (445.0, 45.0, 6, 0.96),
    "H-Alpha": (656.3, 4.0, 2, 0.95),
    "OIII": (500.7, 4.0, 2, 0.94),
    "SII": (672.4, 4.0, 2, 0.93),
}


def theoretical_qe(wavelengths):
    """Gaussian IMX585 QE model in percent (peak ~92% at 550 nm)."""
    return 92.0 * np.exp(-((wavelengths - 550.0) / 200.0) ** 2)


def theoretical_transmission(wavelengths, center, width, power, peak):
    """Analytic filter transmission (percent)."""
    return peak * 100.0 * np.exp(-((wavelengths - center) / width) ** power)


def analyze_spectra(data_dir, config):
    """Load (or synthesize) QE + filter curves and derive summary metrics."""
    a = config.analysis
    wv = np.linspace(a.wavelength_min_nm, a.wavelength_max_nm, a.wavelength_points)
    qe, filters = load_spectra(data_dir)

    if qe is None:
        qe = {"wavelengths": wv, "qe_percentage": theoretical_qe(wv), "synthetic": True}
    if not filters:
        filters = {name: {"wavelengths": wv,
                          "transmission": theoretical_transmission(wv, *params),
                          "synthetic": True}
                   for name, params in FILTER_PROFILES.items()}

    qe_arr = np.asarray(qe["qe_percentage"])
    qe_peak = float(np.max(qe_arr))
    qe_peak_nm = float(np.asarray(qe["wavelengths"])[np.argmax(qe_arr)])

    filter_metrics = {}
    for name, f in filters.items():
        m = _filter_metrics(np.asarray(f["wavelengths"]), np.asarray(f["transmission"]))
        m["synthetic"] = f.get("synthetic", False)
        filter_metrics[name] = m

    return {"qe": qe, "filters": filters, "qe_peak_pct": qe_peak, "qe_peak_nm": qe_peak_nm,
            "qe_synthetic": qe.get("synthetic", False), "filter_metrics": filter_metrics}


def _filter_metrics(wv, trans):
    """Central wavelength, peak transmission, and FWHM of a transmission curve."""
    peak = float(np.max(trans))
    center = float(wv[np.argmax(trans)])
    above = trans >= peak / 2.0
    if above.any():
        idx = np.where(above)[0]
        fwhm = float(wv[idx[-1]] - wv[idx[0]])
    else:
        fwhm = 0.0
    return {"center_nm": center, "peak_pct": peak, "fwhm_nm": fwhm}
