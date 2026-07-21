import argparse
import os
import glob
from astropy.io import fits

# Import the robust engine
from phase3_photometry import UniversalPhotometryEngine

def detect_filter(file_path, default_filter):
    """
    Attempts to read the filter from the FITS header or filename.
    Passes the exact scientific filter (R, V, B) directly to the engine 
    so the Priority Queue can handle the catalog routing!
    """
    try:
        with fits.open(file_path) as hdul:
            header = hdul[0].header
            filt = str(header.get('FILTER', '')).strip().upper()
            
            # If the header is empty, try guessing from the filename
            if not filt or filt == 'NONE':
                base = os.path.basename(file_path).upper()
                if '_R_' in base or '-R-' in base: filt = 'R'
                elif '_V_' in base or '-V-' in base: filt = 'V'
                elif '_B_' in base or '-B-' in base: filt = 'B'
                elif '_I_' in base or '-I-' in base: filt = 'I'
                elif '_G_' in base or '-G-' in base: filt = 'G'
            
            # Map standard observatory filters to strict science bands
            if filt in ['R', 'R-BAND', 'RMAG', 'R_SDSS', 'SDSS-R']: 
                return 'R'
            elif filt in ['I', 'I-BAND', 'IMAG', 'I_SDSS', 'SDSS-I']: 
                return 'I'
            elif filt in ['G', 'G-BAND', 'GMAG', 'G_SDSS', 'SDSS-G']: 
                return 'G'
            elif filt == 'V':
                return 'V'
            elif filt == 'B':
                return 'B'
            else:
                return default_filter
    except Exception:
        # If the file is unreadable or missing headers, fall back to what the user typed
        return default_filter


def main():
    parser = argparse.ArgumentParser(description="Phase 3: Batch Absolute Photometry, Flux Calibration, and Catalog Generation.")
    
    # Required Argument (Accepts a file OR a folder)
    parser.add_argument("input_path", help="Path to a single Master FITS file OR a directory containing Master FITS files.")
    
    # Optional Arguments
    parser.add_argument("--filter", default="R", help="Fallback scientific filter band (e.g., R, G, B, V) if auto-detect fails. Default: R")
    parser.add_argument("--fwhm", type=float, default=3.5, help="Estimated FWHM in pixels. Default: 3.5")
    parser.add_argument("--threshold", type=float, default=5.0, help="Detection threshold in sigma. Default: 5.0")
    parser.add_argument("--outdir", default=None, help="Output directory. Defaults to the directory of the input file(s).")

    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    
    # 1. Determine what files to process (Single file vs. Batch directory)
    files_to_process = []
    
    if os.path.isdir(input_path):
        print(f"Scanning directory for Master FITS files: {input_path}")
        # Find all fits files, safely ignoring ones we've already calibrated or generated
        for f in os.listdir(input_path):
            if f.endswith('.fits') and not f.endswith('_fluxcal.fits') and not f.endswith('_segmap.fits'):
                files_to_process.append(os.path.join(input_path, f))
                
        if not files_to_process:
            print("No valid Master FITS files found in the directory.")
            return
            
        # Sort them alphabetically so they process in a predictable order
        files_to_process.sort() 
        
    elif os.path.isfile(input_path) and input_path.endswith('.fits'):
        # Just a single file was provided
        files_to_process.append(input_path)
    else:
        print("[ERROR] Input path must be a .fits file or a directory containing .fits files.")
        return

    print(f"======================================================")
    print(f"  PHASE 3: BATCH PHOTOMETRY & CATALOG ENGINE")
    print(f"  Found {len(files_to_process)} file(s) to process.")
    print(f"======================================================")

    # 2. Loop through all discovered files
    for i, file_path in enumerate(files_to_process, 1):
        base_name = os.path.basename(file_path).replace('.fits', '')
        input_dir = os.path.dirname(file_path)
        
        # Determine where to save the files
        out_dir = args.outdir if args.outdir else input_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_calibrated_fits = os.path.join(out_dir, f"{base_name}_fluxcal.fits")
        out_catalog = os.path.join(out_dir, f"{base_name}_catalog.csv")
        out_segmap = os.path.join(out_dir, f"{base_name}_segmap.fits")

        # Auto-detect the strict scientific filter for this specific image!
        target_filter = detect_filter(file_path, args.filter)

        print(f"\n[{i}/{len(files_to_process)}] Processing: {base_name}")
        print(f"Target Filter Band: {target_filter} (Auto-detected)")
        
        # VERY IMPORTANT: Re-initialize the Engine for EACH file.
        engine = UniversalPhotometryEngine(fwhm_estimate=args.fwhm, detection_threshold=args.threshold)

        try:
            # Step 1: Calculate Zero Point (Engine will auto-route to APASS/Pan-STARRS)
            engine.calculate_local_zero_point(master_science_fits=file_path, filter_band=target_filter)
            
            # Step 2: Export Calibrated Image (Janskys)
            engine.export_flux_calibrated_image(target_fits=file_path, output_filename=out_calibrated_fits)
            
            # Step 3: Export Catalog and Segmentation Map
            engine.generate_full_catalog(target_fits=file_path, output_csv=out_catalog, output_segmap=out_segmap)
            
            print(f" -> Successfully finished {base_name}")
            
        except Exception as e:
            print(f" -> [ERROR] Failed on {base_name}: {str(e)}")

    # 3. Final Summary
    print("\n--- BATCH PROCESSING COMPLETE ---")
    out_msg = args.outdir if args.outdir else "the respective input directories"
    print(f"All results saved to: {out_msg}")


if __name__ == "__main__":
    main()