"""Characterization orchestrator: data -> PTC/dark/spectral -> results + dashboard."""

import os

import numpy as np

from cassa_camchar.config import load_config
from cassa_camchar.logging_utils import get_logger
from cassa_camchar.io import load_grouped, robust_mean
from cassa_camchar.ptc import analyze_ptc, read_noise_electrons, dynamic_range_stops
from cassa_camchar.dark import temperature_sweep, linearity as dark_linearity
from cassa_camchar.spectral import analyze_spectra
from cassa_camchar.results import CharacterizationResults
from cassa_camchar.report import generate_dashboard


def analyze(data_dir=None, config=None, out_dir=None, logger=None, show=False):
    """Run the full characterization and write results + dashboard.

    Returns the :class:`CharacterizationResults`, or None if there is no data.
    """
    config = config or load_config()
    logger = logger or get_logger("cassa_camchar")
    data_dir = data_dir or config.paths.data_dir
    out_dir = out_dir or config.paths.out_dir
    a = config.analysis
    os.makedirs(out_dir, exist_ok=True)

    logger.info(f"Scanning {data_dir} ...")
    bias, flats, darks = load_grouped(data_dir)
    if not flats:
        logger.error("No flat frames found; cannot run the PTC analysis.")
        return None

    results = CharacterizationResults()
    master_ptc = None

    # --- PTC + gain sweep ---
    for g in sorted(flats.keys()):
        ptc = analyze_ptc(flats[g], config)
        if master_ptc is None and ptc["means"].size:
            master_ptc = ptc
        rn = np.nan
        if g in bias and len(bias[g]) >= 2 and np.isfinite(ptc["gain"]):
            rn = read_noise_electrons(bias[g][0], bias[g][1], ptc["gain"], config)
        dr = dynamic_range_stops(ptc["full_well_e"], rn)
        results.driver_gains.append(g)
        results.system_gain.append(ptc["gain"])
        results.read_noise_e.append(rn)
        results.full_well_e.append(ptc["full_well_e"])
        results.dynamic_range.append(dr)
        logger.info(f"gain {g:>3}: G={ptc['gain']:.3f} e-/ADU | RN={rn:.2f} e- | "
                    f"FWC={ptc['full_well_e']:.0f} e- | DR={dr:.1f} stops")

    # --- Dark current + linearity (reference gain/bias from the lowest driver gain) ---
    ref_gain = next((x for x in results.system_gain if np.isfinite(x)), 1.0)
    ref_bias = robust_mean(bias[sorted(bias)[0]][0], a.roi_fraction, a.sigma_clip) if bias else 0.0
    temps, currents = temperature_sweep(darks, ref_bias, ref_gain, config)
    results.temperatures_c = temps.tolist()
    results.dark_current = currents.tolist()
    lin = dark_linearity(darks, ref_bias, ref_gain, config)
    if lin is not None:
        results.linearity = {"temperature_c": lin["temperature_c"], "slope": lin["slope"],
                             "nonlinearity_pct": lin["nonlinearity_pct"]}
        logger.info(f"dark linearity @ {lin['temperature_c']:.0f}C: "
                    f"{lin['nonlinearity_pct']:.1f}% deviation from linear")

    # --- Spectral (QE + filters) ---
    spectra = analyze_spectra(data_dir, config)
    results.qe_peak_pct = spectra["qe_peak_pct"]
    results.qe_peak_nm = spectra["qe_peak_nm"]
    results.qe_synthetic = spectra["qe_synthetic"]
    results.filter_metrics = spectra["filter_metrics"]
    if spectra["qe_synthetic"]:
        logger.warning("QE curve is synthetic (no lab CSV found) - labelled accordingly.")

    # --- Persist + dashboard ---
    csv_path = os.path.join(out_dir, "characterization_results.csv")
    json_path = os.path.join(out_dir, "characterization_results.json")
    png_path = os.path.join(out_dir, "CASSA_Sensor_Characterization_Report.png")
    results.write_csv(csv_path)
    results.write_json(json_path)
    generate_dashboard(master_ptc, results, spectra, lin, png_path, show=show)
    logger.info(f"Wrote {csv_path}, {json_path}, {png_path}")
    return results
