"""Virtual-camera calibration campaign generator.

Writes synthetic bias/flat/dark FITS frames and QE/filter CSVs on the config's
capture grid. The sensor model below is the ground truth the analysis should
recover, which makes it a validation harness (e.g. system gain at driver gain 0
is 4.0 e-/ADU).
"""

import os
import shutil

import numpy as np
from astropy.io import fits

from cassa_camchar.config import load_config
from cassa_camchar.logging_utils import get_logger
from cassa_camchar.spectral import FILTER_PROFILES, theoretical_qe


# --- Ground-truth sensor model ---------------------------------------------- #
def system_gain(driver_gain):
    """True system gain (e-/ADU): 4.0 at gain 0, decreasing with driver gain."""
    return 4.0 * np.exp(-0.015 * driver_gain)


def read_noise(driver_gain):
    """True read noise (e-): 4.5 at gain 0, decreasing with driver gain."""
    return 3.5 * np.exp(-0.012 * driver_gain) + 1.0


def dark_current(temp_c):
    """True dark current (e-/pixel/s), doubling every ~6 degrees."""
    return 0.002 * (2.0 ** ((temp_c + 20.0) / 6.0))


def _frame(mean_e, rn_e, gain, config, rng):
    size = tuple(config.sensor.image_size)
    shot = rng.normal(mean_e, np.sqrt(max(mean_e, 1e-3)), size)
    read = rng.normal(0, rn_e, size)
    adu = (shot + read) / gain + config.sensor.base_bias_adu
    return np.clip(adu, 0, config.sensor.max_adu).astype(np.uint16)


def _save(path, data, header):
    hdu = fits.PrimaryHDU(data)
    for k, v in header.items():
        hdu.header[k] = v
    hdu.writeto(path, overwrite=True)


def simulate_campaign(out_dir, config=None, seed=None, logger=None):
    """Generate a full synthetic calibration campaign under ``out_dir``."""
    config = config or load_config()
    logger = logger or get_logger("cassa_camchar_sim")
    rng = np.random.default_rng(seed)
    cap = config.capture

    for folder in ("bias", "flat", "dark", "spectra"):
        path = os.path.join(out_dir, folder)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path)

    logger.info("Simulating bias frames...")
    for g in cap.driver_gains:
        rn, sg = read_noise(g), system_gain(g)
        for i in range(2):
            _save(os.path.join(out_dir, "bias", f"bias_gain{g}_{i+1}.fit"),
                  _frame(0, rn, sg, config, rng),
                  {"GAIN": g, "EXPTIME": 0.001, "CCD-TEMP": -10.0, "IMAGETYP": "Bias"})

    logger.info("Simulating flat frames...")
    exposures = np.linspace(cap.exposure_min_s, cap.exposure_max_s, cap.exposure_steps)
    for g in cap.driver_gains:
        rn, sg = read_noise(g), system_gain(g)
        for exp in exposures:
            for i in range(2):
                _save(os.path.join(out_dir, "flat", f"flat_gain{g}_exp{exp:.2f}s_{i+1}.fit"),
                      _frame(cap.flat_electrons_per_s * exp, rn, sg, config, rng),
                      {"GAIN": g, "EXPTIME": round(float(exp), 3), "CCD-TEMP": -10.0, "IMAGETYP": "Flat"})

    logger.info("Simulating dark frames...")
    rn, sg = read_noise(0), system_gain(0)
    for t in cap.temperatures_c:
        _save(os.path.join(out_dir, "dark", f"dark_temp{t}c_exp{int(cap.dark_temp_exposure_s)}s.fit"),
              _frame(dark_current(t) * cap.dark_temp_exposure_s, rn, sg, config, rng),
              {"GAIN": 0, "EXPTIME": float(cap.dark_temp_exposure_s), "CCD-TEMP": float(t), "IMAGETYP": "Dark"})
    dc0 = dark_current(cap.linearity_temp_c)
    for exp in cap.linearity_exposures_s:
        _save(os.path.join(out_dir, "dark", f"dark_lin_temp{int(cap.linearity_temp_c)}c_exp{exp}s.fit"),
              _frame(dc0 * exp, rn, sg, config, rng),
              {"GAIN": 0, "EXPTIME": float(exp), "CCD-TEMP": float(cap.linearity_temp_c), "IMAGETYP": "Dark"})

    logger.info("Simulating spectral CSVs...")
    _write_spectra(os.path.join(out_dir, "spectra"), config, rng)
    logger.info(f"Simulated campaign written to {out_dir}")
    return out_dir


def _write_spectra(spectra_dir, config, rng):
    a = config.analysis
    wv = np.linspace(a.wavelength_min_nm, a.wavelength_max_nm, a.wavelength_points)
    np.savetxt(os.path.join(spectra_dir, "qe_calibration.csv"),
               np.column_stack((wv, theoretical_qe(wv))),
               delimiter=",", header="Wavelength(nm),QE(%)", comments="")
    unfiltered = 50000.0 * np.ones_like(wv)
    for name, (c, w, p, peak) in FILTER_PROFILES.items():
        trans = peak * np.exp(-((wv - c) / w) ** p)
        filtered = np.clip(unfiltered * trans + rng.normal(0, 80, len(wv)), 0, None)
        np.savetxt(os.path.join(spectra_dir, f"filter_{name}.csv"),
                   np.column_stack((wv, unfiltered, filtered)),
                   delimiter=",", header="Wavelength(nm),Unfiltered_Signal,Filtered_Signal", comments="")
