import argparse
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
import os
import glob
from astropy.io import fits

# Catalog APIs
from phase3_photometry import UniversalPhotometryEngine
from astroquery.sdss import SDSS
from astroquery.vizier import Vizier

def detect_filter(csv_path, default_filter):
    """
    Attempts to find the companion FITS file for the CSV, reads the exact FILTER 
    from the FITS header, and maps it. Falls back to filename parsing if FITS is missing.
    """
    base_name_no_ext = os.path.basename(csv_path).replace('_catalog.csv', '')
    input_dir = os.path.dirname(csv_path)
    
    # Look for the corresponding FITS file in the same folder
    possible_fits = [
        os.path.join(input_dir, f"{base_name_no_ext}.fits"),
        os.path.join(input_dir, f"{base_name_no_ext}_fluxcal.fits")
    ]
    
    for pf in possible_fits:
        if os.path.exists(pf):
            try:
                with fits.open(pf) as hdul:
                    header = hdul[0].header
                    filt = str(header.get('FILTER', '')).strip().upper()
                    
                    if filt and filt != 'NONE':
                        # Map standard observatory filters to Pan-STARRS/SDSS catalog bands
                        if filt in ['R', 'R-BAND', 'RMAG', 'R_SDSS', 'SDSS-R']: return 'rmag'
                        elif filt in ['I', 'I-BAND', 'IMAG', 'I_SDSS', 'SDSS-I']: return 'imag'
                        elif filt in ['G', 'G-BAND', 'GMAG', 'G_SDSS', 'SDSS-G']: return 'gmag'
                        elif filt in ['Z', 'Z-BAND', 'ZMAG', 'SDSS-Z']: return 'zmag'
                        elif filt in ['Y', 'Y-BAND', 'YMAG']: return 'ymag'
                        elif filt == 'V': return 'Vmag'
                        elif filt == 'B': return 'Bmag'
            except Exception:
                pass # If FITS reading fails, silently fall back to filename checking

    # Fallback: Guess from the CSV filename if no FITS file exists
    base_name = os.path.basename(csv_path).upper()
    if '_R_' in base_name or '-R-' in base_name: return 'rmag'
    elif '_V_' in base_name or '-V-' in base_name: return 'Vmag'
    elif '_B_' in base_name or '-B-' in base_name: return 'Bmag'
    elif '_I_' in base_name or '-I-' in base_name: return 'imag'
    elif '_G_' in base_name or '-G-' in base_name: return 'gmag'
    
    return default_filter


def cross_match_and_report(catalog_name, cat_coords, cat_mags, bright_stars, my_coords, matched_filter):
    """
    Cross-matches our local stars against a global catalog, calculates the photometric
    error (delta), and prints a formatted scientific report.
    """
    idx, d2d, _ = match_coordinates_sky(my_coords, cat_coords)
    
    # Only accept matches within 2.0 arcseconds to ensure we are looking at the exact same star
    match_mask = d2d < 2.0 * u.arcsec
    
    print(f"\n--- {catalog_name.upper()} VERIFICATION REPORT (Filter: {matched_filter}) ---")
    print(f"{'Obj ID':<8} | {'RA (deg)':<10} | {'DEC (deg)':<10} | {'Your Mag':<10} | {'True Mag':<10} | {'Error (Delta)':<12}")
    print("-" * 75)
    
    errors = []
    for i in range(len(bright_stars)):
        if match_mask[i]:
            ps_idx = idx[i]
            true_mag = cat_mags[ps_idx]
            calc_mag = bright_stars['Absolute_Mag'][i]
            
            # Skip masked, negative, or missing data from the catalog
            if np.ma.is_masked(true_mag) or np.isnan(true_mag) or true_mag < 0: 
                continue
            
            # Calculate the mathematical offset
            delta = calc_mag - true_mag
            errors.append(delta)
            
            print(f"{bright_stars['ID'][i]:<8} | {bright_stars['RA_deg'][i]:.5f} | {bright_stars['Dec_deg'][i]:.5f} | "
                  f"{calc_mag:.3f}      | {true_mag:.3f}      | {delta:+.3f} mag")
                  
    if len(errors) > 0:
        mean_error = np.mean(errors)
        std_error = np.std(errors)
        print("-" * 75)
        print(f"Mean Calibration Error (Offset): {mean_error:+.4f} mag")
        print(f"Standard Deviation (Scatter)   : {std_error:.4f} mag")
        
        # Intelligent feedback based on the catalog
        if abs(mean_error) < 0.1:
            print(f"✅ EXCELLENT! Matches {catalog_name} perfectly.")
        else:
            print(f"⚠️ NOTE: Noticeable offset against {catalog_name}.")
            print(f"   This is likely a 'Color Term' difference between your physical glass filter")
            print(f"   and the {catalog_name} passband. To force-match this specific catalog,")
            print(f"   you would need to add {-mean_error:.4f} to your Zero Point.")
    else:
        print(f"No valid matches found between your stars and {catalog_name}.")


def verify_calibration(csv_path, filter_band="rmag"):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Cannot find catalog file: {csv_path}")

    print(f"Analyzing Catalog: {os.path.basename(csv_path)}")
    print(f"Target Filter Band: {filter_band}")
    
    # 1. Load the generated catalog
    cat = Table.read(csv_path, format='csv')
    
    # 2. Filter for STARS (Round objects)
    # An ellipticity of < 0.15 ensures we drop galaxies, cosmic rays, and hot pixels
    if 'Ellipticity' not in cat.colnames:
        print("[ERROR] 'Ellipticity' column missing. Cannot separate stars from galaxies.")
        return
        
    mask_stars = cat['Ellipticity'] < 0.15
    stars = cat[mask_stars]
    
    if len(stars) == 0:
        print("[ERROR] No round stars found in the catalog based on the ellipticity threshold.")
        return

    # 3. Sort by brightness and take the top 25
    # We only use the brightest stars because they have the highest Signal-to-Noise ratio
    stars.sort('Absolute_Mag')
    bright_stars = stars[:25]
    
    print(f" -> Found {len(stars)} round objects (stars). Selecting the {len(bright_stars)} brightest for verification...\n")
    
    # 4. Get the center coordinate and dynamic bounding box
    mean_ra = np.mean(bright_stars['RA_deg'])
    mean_dec = np.mean(bright_stars['Dec_deg'])
    center_coord = SkyCoord(ra=mean_ra, dec=mean_dec, unit=(u.deg, u.deg))
    my_coords = SkyCoord(ra=bright_stars['RA_deg'], dec=bright_stars['Dec_deg'], unit=(u.deg, u.deg))
    
    ra_spread = np.max(bright_stars['RA_deg']) - np.min(bright_stars['RA_deg'])
    dec_spread = np.max(bright_stars['Dec_deg']) - np.min(bright_stars['Dec_deg'])
    radius = max(ra_spread, dec_spread) * u.deg / 2.0 + (1 * u.arcmin)
    
    # =================================================================
    # CATALOG 1: Pan-STARRS (via NASA MAST)
    # =================================================================
    try:
        if filter_band in ['Vmag', 'Bmag']:
            print(f">>> Skipping Pan-STARRS DR2 - Filter '{filter_band}' is not supported by Pan-STARRS.")
        else:
            print(f">>> Querying Pan-STARRS DR2 (NASA MAST) for radius {radius.to(u.arcmin):.2f}...")
            engine = UniversalPhotometryEngine()
            panstarrs_cat = engine._query_panstarrs(center_coord, radius)
            ps_coords = SkyCoord(ra=panstarrs_cat['RAJ2000'], dec=panstarrs_cat['DEJ2000'], unit=(u.deg, u.deg))
            ps_mags = panstarrs_cat[filter_band]
            cross_match_and_report("Pan-STARRS", ps_coords, ps_mags, bright_stars, my_coords, filter_band)
    except Exception as e:
        print(f"[WARNING] Pan-STARRS Query Failed: {e}")

    # =================================================================
    # CATALOG 2: SDSS - Sloan Digital Sky Survey
    # =================================================================
    try:
        if filter_band in ['Vmag', 'Bmag']:
            print(f"\n>>> Skipping SDSS DR12 - Filter '{filter_band}' is not supported by SDSS.")
        else:
            print("\n>>> Querying SDSS DR12 (via SQL to bypass radius limits)...")
            # Convert 'rmag' -> 'r' for SDSS formatting
            sdss_filter = filter_band[0] if filter_band.endswith('mag') else filter_band 
            
            # Build a bounding box to strictly encapsulate our 25 stars
            ra_min, ra_max = np.min(bright_stars['RA_deg']) - 0.02, np.max(bright_stars['RA_deg']) + 0.02
            dec_min, dec_max = np.min(bright_stars['Dec_deg']) - 0.02, np.max(bright_stars['Dec_deg']) + 0.02
            
            # Send a raw SQL query to SDSS to bypass the arbitrary 3.0 arcmin Python limit (type=6 means 'STAR')
            sql_query = f"SELECT ra, dec, {sdss_filter} FROM PhotoObj WHERE ra BETWEEN {ra_min} AND {ra_max} AND dec BETWEEN {dec_min} AND {dec_max} AND type=6"
            sdss_res = SDSS.query_sql(sql_query)
            
            if sdss_res is not None and len(sdss_res) > 0:
                sdss_coords = SkyCoord(ra=sdss_res['ra'], dec=sdss_res['dec'], unit=(u.deg, u.deg))
                sdss_mags = sdss_res[sdss_filter]
                cross_match_and_report("SDSS", sdss_coords, sdss_mags, bright_stars, my_coords, sdss_filter)
            else:
                print(" -> Target field is outside the SDSS sky footprint or no stars matched.")
    except Exception as e:
        print(f"[WARNING] SDSS Query Failed: {e}")

    # =================================================================
    # CATALOG 3: APASS - AAVSO Photometric All-Sky Survey
    # =================================================================
    try:
        print("\n>>> Querying APASS DR9 (VizieR)...")
        
        # We use an Asian-focused robust mirror hopper to bypass ISP blocks
        apass_res = None
        mirrors = [
            'vizier.iucaa.in',         # India (Pune - Best chance for South Asia)
            'vizier.nao.ac.jp',        # Japan (Highly reliable fallback)
            'vizier.china-vo.org',     # China
            'vizier.saao.ac.za',       # South Africa
            'vizier.cds.unistra.fr',   # France
            'vizier.cfa.harvard.edu',  # USA
            'vizier.ast.cam.ac.uk',    # UK 
            'vizier.hia.nrc.ca'        # Canada
        ]
        
        for mirror in mirrors:
            print(f" -> Attempting APASS connection via mirror: {mirror}")
            # Request all columns (*) instead of risking a specific name
            v = Vizier(columns=['RAJ2000', 'DEJ2000', '*'], row_limit=500)
            v.VIZIER_SERVER = mirror
            v.TIMEOUT = 15
            try:
                res = v.query_region(center_coord, radius=radius, catalog='II/336/apass9')
                if len(res) > 0:
                    apass_res = res[0]
                    print(f" -> Success! Connected to {mirror}")
                    break # Success! Break out of the mirror loop
            except Exception as e:
                print(f" -> [WARNING] Mirror {mirror} failed or blocked.")
        
        if apass_res is not None:
            # Bulletproof string stripping to exactly match the Engine's fuzzy logic
            filter_prefix = filter_band[0].lower() # 'r', 'v', 'b'
            
            target_cols = [
                c for c in apass_res.colnames 
                if c.replace("'", "").replace("_", "").lower().startswith(filter_prefix) 
                and 'mag' in c.lower()
            ]
            
            if len(target_cols) > 0:
                best_col = target_cols[0]
                print(f" -> Auto-detected APASS magnitude column: '{best_col}'")
                apass_coords = SkyCoord(ra=apass_res['RAJ2000'], dec=apass_res['DEJ2000'], unit=(u.deg, u.deg))
                apass_mags = apass_res[best_col]
                cross_match_and_report("APASS", apass_coords, apass_mags, bright_stars, my_coords, best_col)
            else:
                print(f" -> Connected, but could not find a column matching filter '{filter_band}'. Available columns: {apass_res.colnames}")
        else:
            print(" -> No APASS stars found or all global VizieR mirrors were blocked.")
    except Exception as e:
        print(f"[WARNING] APASS process failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Verify Photometric Calibration against Multiple Global Catalogs")
    parser.add_argument("input_path", help="Path to a single _catalog.csv file OR a directory containing _catalog.csv files.")
    parser.add_argument("--filter", default="rmag", help="Fallback filter band if auto-detect fails (e.g., rmag). Default: rmag")
    args = parser.parse_args()
    
    input_path = os.path.abspath(args.input_path)
    csv_files = []
    
    # 1. Determine files to process
    if os.path.isdir(input_path):
        print(f"Scanning directory for Catalog CSV files: {input_path}")
        for f in os.listdir(input_path):
            if f.endswith('_catalog.csv'):
                csv_files.append(os.path.join(input_path, f))
        csv_files.sort()
        if not csv_files:
            print("No valid _catalog.csv files found in the directory.")
            return
    elif os.path.isfile(input_path) and input_path.endswith('.csv'):
        csv_files.append(input_path)
    else:
        print("[ERROR] Input path must be a .csv file or a directory containing _catalog.csv files.")
        return

    print(f"======================================================")
    print(f"  MULTI-CATALOG CALIBRATION VERIFICATION TOOL (BATCH)")
    print(f"  Found {len(csv_files)} file(s) to verify.")
    print(f"======================================================")

    # 2. Iterate through all found CSVs
    for i, csv_file in enumerate(csv_files, 1):
        
        # Use the robust FITS header check
        target_filter = detect_filter(csv_file, args.filter)
        
        print(f"\n\n{'='*75}")
        print(f"[{i}/{len(csv_files)}] VERIFYING FILE")
        print(f"{'='*75}")
        
        verify_calibration(csv_file, target_filter)
        
    print("\n--- BATCH VERIFICATION COMPLETE ---")

if __name__ == "__main__":
    main()