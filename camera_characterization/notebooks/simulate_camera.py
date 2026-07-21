import numpy as np
from astropy.io import fits
import os
import shutil

# ==========================================
# 1. VIRTUAL CAMERA SPECIFICATIONS
# ==========================================
IMAGE_SIZE = (500, 500)  
BASE_BIAS_ADU = 500      
MAX_ADU = 65535          

def get_system_gain(driver_gain): return 4.0 * np.exp(-0.015 * driver_gain)
def get_read_noise(driver_gain): return 3.5 * np.exp(-0.012 * driver_gain) + 1.0
def get_dark_current(temp_c): return 0.002 * (2.0 ** ((temp_c + 20.0) / 6.0))

def generate_frame(mean_signal_electrons, read_noise_electrons, system_gain):
    shot_noise = np.random.normal(mean_signal_electrons, np.sqrt(max(mean_signal_electrons, 0.001)), IMAGE_SIZE)
    read_noise = np.random.normal(0, read_noise_electrons, IMAGE_SIZE)
    frame_adu = ((shot_noise + read_noise) / system_gain) + BASE_BIAS_ADU
    return np.clip(frame_adu, 0, MAX_ADU).astype(np.uint16)

def save_fits(filename, data, header_info):
    hdu = fits.PrimaryHDU(data)
    for k, v in header_info.items(): hdu.header[k] = v
    hdu.writeto(filename, overwrite=True)

# ==========================================
# 2. DATA GENERATION PIPELINE
# ==========================================
def simulate_camera_campaign(base_dir="./"):
    print("Initializing Virtual Camera Calibration Campaign...")
    
    folders = ['bias', 'flat', 'dark', 'spectra']
    for folder in folders:
        path = os.path.join(base_dir, folder)
        if os.path.exists(path): shutil.rmtree(path)
        os.makedirs(path)
        
    driver_gains = [0, 50, 100, 150, 200]
    
    # ---------------------------------------------------------
    # A. BIAS
    # ---------------------------------------------------------
    print("Shooting Bias Frames...")
    for gain in driver_gains:
        rn, sys_gain = get_read_noise(gain), get_system_gain(gain)
        for i in range(2): 
            save_fits(os.path.join(base_dir, 'bias', f'bias_gain{gain}_{i+1}.fit'), 
                      generate_frame(0, rn, sys_gain), 
                      {'GAIN': gain, 'EXPTIME': 0.001, 'CCD-TEMP': -10.0, 'IMAGETYP': 'Bias'})

    # ---------------------------------------------------------
    # B. FLAT
    # ---------------------------------------------------------
    print("Shooting Flat Frames...")
    for gain in driver_gains:
        rn, sys_gain = get_read_noise(gain), get_system_gain(gain)
        for exp in np.linspace(0.1, 12.0, 20):
            for i in range(2): 
                save_fits(os.path.join(base_dir, 'flat', f'flat_gain{gain}_exp{exp:.2f}s_{i+1}.fit'), 
                          generate_frame(5000.0 * exp, rn, sys_gain), 
                          {'GAIN': gain, 'EXPTIME': round(exp, 3), 'CCD-TEMP': -10.0, 'IMAGETYP': 'Flat'})

    # ---------------------------------------------------------
    # C. DARK (Temp Sweep)
    # ---------------------------------------------------------
    print("Shooting Dark Frames (Thermal Sweep)...")
    rn, sys_gain = get_read_noise(0), get_system_gain(0)
    for temp in [-20, -15, -10, -5, 0, 5, 10, 15, 20]:
        save_fits(os.path.join(base_dir, 'dark', f'dark_temp{temp}c_exp300s.fit'), 
                  generate_frame(get_dark_current(temp) * 300.0, rn, sys_gain), 
                  {'GAIN': 0, 'EXPTIME': 300.0, 'CCD-TEMP': float(temp), 'IMAGETYP': 'Dark'})

    # ---------------------------------------------------------
    # D. DARK (Time Sweep)
    # ---------------------------------------------------------
    print("Shooting Dark Linearity Frames...")
    dc_rate = get_dark_current(0.0)
    for exp in [10.0, 30.0, 60.0, 120.0, 300.0, 600.0]:
        save_fits(os.path.join(base_dir, 'dark', f'dark_temp0c_exp{exp}s.fit'), 
                  generate_frame(dc_rate * exp, rn, sys_gain), 
                  {'GAIN': 0, 'EXPTIME': float(exp), 'CCD-TEMP': 0.0, 'IMAGETYP': 'Dark'})

    # ---------------------------------------------------------
    # E. SPECTRAL CSVs (LRGBSHO Filters)
    # ---------------------------------------------------------
    print("Exporting Spectral CSVs (QE & LRGBSHO Filters)...")
    spectra_dir = os.path.join(base_dir, 'spectra')
    wavelengths = np.linspace(350, 1000, 500) # High res array for smooth curves
    
    # 1. Quantum Efficiency
    np.savetxt(os.path.join(spectra_dir, 'qe_calibration.csv'), 
               np.column_stack((wavelengths, 92 * np.exp(-((wavelengths - 550) / 200) ** 2))), 
               delimiter=',', header='Wavelength(nm),QE(%)', comments='')
               
    # 2. LRGBSHO Filter Profiles
    # Broadband uses Power=6 for a "flat-top" dichroic look. Narrowband uses Power=2 for a Gaussian peak.
    filter_profiles = {
        "Luminance": {"center": 550.0, "width": 150.0, "power": 6, "peak": 0.98},
        "Red":       {"center": 640.0, "width": 55.0,  "power": 6, "peak": 0.97},
        "Green":     {"center": 535.0, "width": 45.0,  "power": 6, "peak": 0.96},
        "Blue":      {"center": 445.0, "width": 45.0,  "power": 6, "peak": 0.96},
        "H-Alpha":   {"center": 656.3, "width": 4.0,   "power": 2, "peak": 0.95},
        "OIII":      {"center": 500.7, "width": 4.0,   "power": 2, "peak": 0.94},
        "SII":       {"center": 672.4, "width": 4.0,   "power": 2, "peak": 0.93}
    }
    
    unfiltered_signal = 50000 * np.ones_like(wavelengths) # Flat baseline light source
    
    for name, params in filter_profiles.items():
        c = params["center"]
        w = params["width"]
        p = params["power"]
        peak = params["peak"]
        
        # Mathematical model of the transmission curve
        ideal_trans = peak * np.exp(-((wavelengths - c) / w) ** p) 
        
        # Apply the transmission model to the light source and add realistic lab noise
        filtered_signal = np.clip(unfiltered_signal * ideal_trans + np.random.normal(0, 80, len(wavelengths)), 0, None)
        
        np.savetxt(os.path.join(spectra_dir, f'filter_{name}.csv'), 
                   np.column_stack((wavelengths, unfiltered_signal, filtered_signal)), 
                   delimiter=',', header='Wavelength(nm),Unfiltered_Signal,Filtered_Signal', comments='')

    print("\nSimulation Complete!")

if __name__ == "__main__":
    TARGET_DIRECTORY = "../data" 
    simulate_camera_campaign(base_dir=TARGET_DIRECTORY)