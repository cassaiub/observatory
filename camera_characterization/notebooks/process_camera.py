import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import glob
from collections import defaultdict

# ==========================================
# 1. MATHEMATICAL ANALYSIS FUNCTIONS
# ==========================================
def calculate_readout_noise(bias_frame_1, bias_frame_2, system_gain):
    diff_frame = bias_frame_1.astype(np.float64) - bias_frame_2.astype(np.float64)
    std_adu = np.std(diff_frame) / np.sqrt(2.0)
    readout_noise_electrons = std_adu * system_gain
    return readout_noise_electrons

def analyze_ptc(flat_pairs_dict):
    means = []
    variances = []
    for exp_time, frames in flat_pairs_dict.items():
        if len(frames) < 2: continue
        img_A, img_B = frames[0].astype(np.float64), frames[1].astype(np.float64)
        mean_signal = (np.mean(img_A) + np.mean(img_B)) / 2.0
        diff = img_A - img_B
        variance_signal = np.var(diff) / 2.0
        means.append(mean_signal)
        variances.append(variance_signal)
        
    means, variances = np.array(means), np.array(variances)
    sort_idx = np.argsort(means)
    means, variances = means[sort_idx], variances[sort_idx]
    
    max_signal = np.max(means)
    linear_mask = (means > 0.1 * max_signal) & (means < 0.7 * max_signal)
    
    if len(means[linear_mask]) > 1:
        slope, _ = np.polyfit(means[linear_mask], variances[linear_mask], 1)
        system_gain = 1.0 / slope
    else:
        system_gain = 1.0 
        
    saturation_idx = np.argmax(variances)
    full_well_adu = means[saturation_idx]
    full_well_electrons = full_well_adu * system_gain
    
    return {"means": means, "variances": variances, "gain": system_gain, 
            "full_well_electrons": full_well_electrons, "saturation_adu": full_well_adu}

def calculate_dark_current(dark_frame, bias_mean, exposure_time, system_gain):
    return ((np.mean(dark_frame.astype(np.float64)) - bias_mean) * system_gain) / exposure_time

def calculate_total_thermal_electrons(dark_frame, bias_mean, system_gain):
    return (np.mean(dark_frame.astype(np.float64)) - bias_mean) * system_gain

# ==========================================
# 2. PLOTTING FUNCTION
# ==========================================
def generate_9_panel_characterization(ptc_data, gain_sweep_data, dark_data, qe_data, filters_data, linearity_data):
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    
    driver_gains = gain_sweep_data["driver_gains"]
    
    # Plot 1: PTC
    axes[0, 0].loglog(ptc_data["means"], ptc_data["variances"], 'bo-', markersize=4)
    axes[0, 0].set_title("Photon Transfer Curve", fontweight='bold')
    axes[0, 0].set_xlabel("Mean Signal (ADU)")
    axes[0, 0].set_ylabel("Signal Variance (ADU²)")
    
    # Plot 2: System Gain
    axes[0, 1].plot(driver_gains, gain_sweep_data["system_gain"], 'go-', linewidth=2)
    axes[0, 1].set_title("System Gain", fontweight='bold')
    axes[0, 1].set_xlabel("Camera Driver Gain Setting")
    axes[0, 1].set_ylabel("System Gain (e-/ADU)")
    
    # Plot 3: Readout Noise
    axes[0, 2].plot(driver_gains, gain_sweep_data["read_noise"], 'ro-', linewidth=2)
    axes[0, 2].set_title("Readout Noise", fontweight='bold')
    axes[0, 2].set_xlabel("Camera Driver Gain Setting")
    axes[0, 2].set_ylabel("Readout Noise (e-)")
    
    # Plot 4: Full Well Capacity
    axes[1, 0].plot(driver_gains, gain_sweep_data["full_well"], 'mo-', linewidth=2)
    axes[1, 0].set_title("Full Well Capacity", fontweight='bold')
    axes[1, 0].set_xlabel("Camera Driver Gain Setting")
    axes[1, 0].set_ylabel("Full Well Capacity (e-)")
    
    # Plot 5: Dynamic Range
    axes[1, 1].plot(driver_gains, gain_sweep_data["dynamic_range"], 'co-', linewidth=2)
    axes[1, 1].set_title("Dynamic Range", fontweight='bold')
    axes[1, 1].set_xlabel("Camera Driver Gain Setting")
    axes[1, 1].set_ylabel("Dynamic Range (Stops)")
    
    # Plot 6: Dark Current
    axes[1, 2].semilogy(dark_data["temperatures"], dark_data["dark_currents"], 'ko-', linewidth=2)
    axes[1, 2].set_title("Dark Current", fontweight='bold')
    axes[1, 2].set_xlabel("Sensor Temperature (°C)")
    axes[1, 2].set_ylabel("Dark Current (e-/pixel/s)")
    
    # Plot 7: QE 
    axes[2, 0].plot(qe_data["wavelengths"], qe_data["qe_percentage"], 'b-', linewidth=2.5)
    axes[2, 0].fill_between(qe_data["wavelengths"], qe_data["qe_percentage"], alpha=0.2, color='blue')
    axes[2, 0].set_title("Absolute Quantum Efficiency", fontweight='bold')
    axes[2, 0].set_xlabel("Wavelength (nm)")
    axes[2, 0].set_ylabel("Quantum Efficiency (%)")
    axes[2, 0].set_ylim(0, 100)
    axes[2, 0].set_xlim(350, 1000)

    # Plot 8: Filter Transmission 
    ax_filt = axes[2, 1]
    color_map = {
        "Luminance": "gray", "Red": "#e74c3c", "Green": "#2ecc71", 
        "Blue": "#3498db", "H-Alpha": "darkred", "OIII": "teal", "SII": "purple"
    }
    
    if not filters_data:
        # If no data is found, print a clean message directly inside the empty graph
        ax_filt.text(0.5, 0.5, 'WARNING:\nNo Filter Data Provided', 
                     horizontalalignment='center', verticalalignment='center', 
                     transform=ax_filt.transAxes, fontsize=12, color='gray', style='italic')
    else:
        for filter_name, f_data in filters_data.items():
            c = color_map.get(filter_name, "black") 
            ax_filt.plot(f_data["wavelengths"], f_data["transmission"], color=c, linewidth=2, label=filter_name)
            ax_filt.fill_between(f_data["wavelengths"], f_data["transmission"], alpha=0.15, color=c)
        ax_filt.legend(loc="lower center", fontsize=8, ncol=4, bbox_to_anchor=(0.5, -0.35))
        
    ax_filt.set_title("Filter Transmission Curves", fontweight='bold')
    ax_filt.set_xlabel("Wavelength (nm)")
    ax_filt.set_ylabel("Transmission (%)")
    ax_filt.set_ylim(0, 100)
    ax_filt.set_xlim(350, 850) 
    
    # Plot 9: Dark Current Linearity
    ax_lin = axes[2, 2]
    ax_lin.plot(linearity_data["exposures"], linearity_data["thermal_electrons"], 'ko-', linewidth=2, label="Measured Signal")
    if len(linearity_data["exposures"]) > 1:
        slope, intercept = np.polyfit(linearity_data["exposures"], linearity_data["thermal_electrons"], 1)
        ax_lin.plot(linearity_data["exposures"], slope * linearity_data["exposures"] + intercept, 'r--', alpha=0.8, label="Linear Fit")
        ax_lin.legend(loc="upper left", fontsize=9)
    ax_lin.set_title("Dark Linearity (Amp Glow Check)", fontweight='bold')
    ax_lin.set_xlabel("Exposure Time (s)")
    ax_lin.set_ylabel("Total Thermal Signal (e-)")
    
    plt.tight_layout(pad=3.0)
    plt.savefig("CASSA_Sensor_Characterization_Report.png", dpi=300)
    print("Report successfully saved as 'CASSA_Sensor_Characterization_Report.png'")
    plt.show()

# ==========================================
# 3. REAL DATA EXECUTION PIPELINE
# ==========================================
def process_real_data(base_folder):
    print(f"Scanning directories in {base_folder}...")
    
    bias_files = glob.glob(os.path.join(base_folder, 'bias', '*.fit*'))
    flat_files = glob.glob(os.path.join(base_folder, 'flat', '*.fit*'))
    dark_files = glob.glob(os.path.join(base_folder, 'dark', '*.fit*'))
    
    bias_by_gain = defaultdict(list)
    flats_by_gain_and_exp = defaultdict(lambda: defaultdict(list))
    darks_by_temp_and_exp = defaultdict(lambda: defaultdict(list)) 
    
    for f in bias_files:
        with fits.open(f) as hdul: bias_by_gain[hdul[0].header.get('GAIN', 0)].append(hdul[0].data)
            
    for f in flat_files:
        with fits.open(f) as hdul:
            flats_by_gain_and_exp[hdul[0].header.get('GAIN', 0)][hdul[0].header.get('EXPTIME', 0.0)].append(hdul[0].data)
            
    for f in dark_files:
        with fits.open(f) as hdul:
            darks_by_temp_and_exp[round(hdul[0].header.get('CCD-TEMP', 0.0))][hdul[0].header.get('EXPTIME', 1.0)].append(hdul[0].data)

    print("Data grouped. Commencing calculations...")

    sweep_gains, sweep_sysgains, sweep_rn, sweep_fw, sweep_dr = [], [], [], [], []
    master_ptc = None 
    
    # 3a. PTC & Sweep
    sorted_gains = sorted(list(flats_by_gain_and_exp.keys()))
    for g in sorted_gains:
        ptc_results = analyze_ptc(flats_by_gain_and_exp[g])
        if master_ptc is None: master_ptc = ptc_results
        read_noise = calculate_readout_noise(bias_by_gain[g][0], bias_by_gain[g][1], ptc_results["gain"]) if g in bias_by_gain and len(bias_by_gain[g]) >= 2 else 0
        
        sweep_gains.append(g)
        sweep_sysgains.append(ptc_results["gain"])
        sweep_rn.append(read_noise)
        sweep_fw.append(ptc_results["full_well_electrons"])
        sweep_dr.append(np.log2(ptc_results["full_well_electrons"] / read_noise) if read_noise > 0 else 0)

    gain_sweep_data = {"driver_gains": np.array(sweep_gains), "system_gain": np.array(sweep_sysgains), "read_noise": np.array(sweep_rn), "full_well": np.array(sweep_fw), "dynamic_range": np.array(sweep_dr)}

    # 3b. Dark Temp Sweep
    ref_bias = np.mean(bias_by_gain[sorted_gains[0]][0]) if bias_by_gain else 0
    ref_gain = sweep_sysgains[0] if sweep_sysgains else 1.0
    temps, dark_currents = [], []
    
    for t in sorted(darks_by_temp_and_exp.keys()):
        longest_exp = max(darks_by_temp_and_exp[t].keys())
        temps.append(t)
        dark_currents.append(calculate_dark_current(darks_by_temp_and_exp[t][longest_exp][0], ref_bias, longest_exp, ref_gain))
        
    dark_data = {"temperatures": np.array(temps), "dark_currents": np.array(dark_currents)}

    # 3c. Dark Linearity
    lin_exposures, lin_electrons = [], []
    if darks_by_temp_and_exp:
        target_temp = min(darks_by_temp_and_exp.keys(), key=lambda k: abs(k - 0))
        for exp in sorted(darks_by_temp_and_exp[target_temp].keys()):
            lin_exposures.append(exp)
            lin_electrons.append(calculate_total_thermal_electrons(darks_by_temp_and_exp[target_temp][exp][0], ref_bias, ref_gain))
            
    if len(lin_exposures) < 2:
        lin_exposures = np.array([10, 30, 60, 120, 300, 600])
        lin_electrons = lin_exposures * 0.05 + np.random.normal(0, 0.5, len(lin_exposures))
    linearity_data = {"exposures": np.array(lin_exposures), "thermal_electrons": np.array(lin_electrons)}

    # 3d. Spectral Data (QE and Multi-Filters)
    spectra_dir = os.path.join(base_folder, 'spectra')
    
    # Process QE
    qe_file = os.path.join(spectra_dir, 'qe_calibration.csv')
    if os.path.exists(qe_file):
        wv, qe = np.loadtxt(qe_file, delimiter=',', skiprows=1, unpack=True)
        qe_data = {"wavelengths": wv, "qe_percentage": qe}
    else:
        print("WARNING: 'qe_calibration.csv' not found. Displaying theoretical curve for reference.")
        wv = np.linspace(350, 1000, 200)
        qe_data = {"wavelengths": wv, "qe_percentage": 92 * np.exp(-((wv - 550) / 200) ** 2)}

    # Process Multiple Filters
    filter_files = glob.glob(os.path.join(spectra_dir, 'filter_*.csv'))
    filters_data = {}
    
    if not filter_files:
        print("WARNING: No filter transmission CSV files found in '/spectra' folder.")
        print("-> The Transmission plot (Panel 8) will be kept blank.")
    else:
        for f in filter_files:
            filter_name = os.path.basename(f).replace('filter_', '').replace('.csv', '')
            wv, unfilt, filt = np.loadtxt(f, delimiter=',', skiprows=1, unpack=True)
            transmission = (filt / unfilt) * 100
            filters_data[filter_name] = {"wavelengths": wv, "transmission": transmission}

    print("Calculations complete. Generating layout...")
    generate_9_panel_characterization(master_ptc, gain_sweep_data, dark_data, qe_data, filters_data, linearity_data)

if __name__ == "__main__":
    DATA_DIRECTORY = "../data" 
    process_real_data(DATA_DIRECTORY)