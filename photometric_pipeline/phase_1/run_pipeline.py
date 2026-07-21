# run_pipeline.py
import os
import glob
import argparse
import numpy as np
from astropy.io import fits
import ccdproc
from astropy.nddata import CCDData
from tqdm import tqdm

# Import our universal network profile
from instruments import ITelescopeNetworkProfile
from data_models import load_standardized_ccds
from processor import UniversalProcessor

def build_master_bias(file_list, instrument, desc="Master Bias"):
    """Builds a Master Bias by median combining raw bias frames."""
    if not file_list: 
        return None
    ccd_list = []
    for f in tqdm(file_list, desc=desc, unit="frame"):
        ccd = load_standardized_ccds(f, instrument)[0].ccd
        ccd_list.append(ccd)
        
    tqdm.write(f" -> Calculating median combine for {len(file_list)} Bias frames...")
    return ccdproc.Combiner(ccd_list).median_combine()

def build_master_dark(file_list, instrument, master_bias, desc="Master Dark"):
    """Builds a Master Dark. Subtracts Master Bias from each frame first."""
    if not file_list: 
        return None
    ccd_list = []
    for f in tqdm(file_list, desc=desc, unit="frame"):
        ccd = load_standardized_ccds(f, instrument)[0].ccd
        if master_bias is not None:
            # Prevent double-bias subtraction by calibrating the darks now
            ccd = ccdproc.subtract_bias(ccd, master_bias)
        ccd_list.append(ccd)
        
    tqdm.write(f" -> Calculating median combine for {len(file_list)} Dark frames...")
    return ccdproc.Combiner(ccd_list).median_combine()

def build_master_flat(file_list, instrument, master_bias, desc="Master Flat"):
    """Builds a Master Flat. Subtracts Master Bias, combines, and NORMALIZES to 1.0."""
    if not file_list: 
        return None
    ccd_list = []
    for f in tqdm(file_list, desc=desc, unit="frame"):
        ccd = load_standardized_ccds(f, instrument)[0].ccd
        if master_bias is not None:
            # Calibrate the raw flat to reveal the true optical shadows
            ccd = ccdproc.subtract_bias(ccd, master_bias)
        ccd_list.append(ccd)
        
    tqdm.write(f" -> Calculating median combine for {len(file_list)} Flat frames...")
    combiner = ccdproc.Combiner(ccd_list)
    m_flat = combiner.median_combine()
    
    # --- CRITICAL FIX: NORMALIZATION ---
    flat_median = np.nanmedian(m_flat.data)
    if flat_median > 0:
        tqdm.write(f" -> Normalizing Master Flat (Dividing by median ADU: {flat_median:.1f})")
        m_flat.data = m_flat.data / flat_median
    # -----------------------------------
    return m_flat

def run(data_dir, output_dir):
    instrument = ITelescopeNetworkProfile()
    print(f"\n--- Initializing Pipeline: {instrument.name} ---")
    print(f"Input Directory:  {data_dir}")
    print(f"Output Directory: {output_dir}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("\nScanning and categorizing files...")
    all_files = glob.glob(os.path.join(data_dir, "*.f*t*")) 
    
    if not all_files:
        print(f"Error: No FITS files found in {data_dir}")
        return

    categorized = {'bias': [], 'dark': [], 'flat': {}, 'science': {}}
    
    for f in tqdm(all_files, desc="Scanning headers", unit="file"):
        try:
            with fits.open(f) as hdul:
                header = hdul[0].header
                img_type = instrument.get_image_type(header)
                
                if img_type in ['flat', 'science']:
                    filt = instrument.get_filter(header)
                    if filt not in categorized[img_type]:
                        categorized[img_type][filt] = []
                    categorized[img_type][filt].append(f)
                else:
                    categorized[img_type].append(f)
        except Exception as e:
            tqdm.write(f"Error reading {os.path.basename(f)}: {e}")

    # 3. Build Master Calibrations
    print("\n--- Building Master Frames ---")
    m_bias = build_master_bias(categorized['bias'], instrument, desc="Master Bias")
    m_dark = build_master_dark(categorized['dark'], instrument, m_bias, desc="Master Dark")
    
    m_flats = {}
    for filt, files in categorized['flat'].items():
        m_flats[filt] = build_master_flat(files, instrument, m_bias, desc=f"Master Flat ({filt})")

    # Proxy flat logic for missing red flat
    if 'RED' not in m_flats and 'LUMINANCE' in m_flats:
        print("\n[Notice] Missing RED flat. Mapping LUMINANCE Master Flat as a proxy for RED.")
        m_flats['RED'] = m_flats['LUMINANCE']

    # 4. Initialize the Universal Processor
    print("\n--- Initializing Processor Engine ---")
    processor = UniversalProcessor(master_bias=m_bias, master_dark=m_dark, master_flats=m_flats)
    
    # 5. Process Science Frames sequentially by filter
    print("\n--- Processing Science Frames ---")
    if not categorized['science']:
        print("No science frames found to process.")
        
    for filt, sci_files in categorized['science'].items():
        print(f"\n>> Processing {filt} filter frames <<")
        
        for sci_file in tqdm(sci_files, desc=f"Calibrating {filt}", unit="img"):
            standard_ccds = load_standardized_ccds(sci_file, instrument)
            calibrated_ccds = processor.process_science_frame(standard_ccds)
            
            for i, cal_ccd in enumerate(calibrated_ccds):
                out_name = f"calibrated_{os.path.basename(sci_file)}"
                out_path = os.path.join(output_dir, out_name)
                cal_ccd.write(out_path, overwrite=True)
            
    print("\n✅ Pipeline execution complete. Your science-ready data is in the output folder.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Phase 1 Instrument Signature Removal Pipeline.")
    parser.add_argument("-i", "--input", required=True, help="Path to the directory containing raw FITS files.")
    parser.add_argument("-o", "--output", required=True, help="Path to the directory where calibrated files will be saved.")
    
    args = parser.parse_args()
    run(args.input, args.output)