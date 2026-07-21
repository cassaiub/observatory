import numpy as np
import time
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from astroquery.vizier import Vizier
from astroquery.sdss import SDSS
from astropy.table import Table
import urllib.request
import json

# Photutils imports for Point Sources (Zero Point Calculation)
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry

# Photutils imports for Extended Objects (Galaxies/Full Catalog)
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
from astropy.convolution import Gaussian2DKernel, convolve

import ois 

class UniversalPhotometryEngine:
    def __init__(self, fwhm_estimate=3.5, detection_threshold=5.0):
        self.fwhm = fwhm_estimate
        self.threshold = detection_threshold
        self.zero_point = None

    def _query_apass(self, center_coord, radius):
        """
        PRIORITY 1 CATALOG: APASS DR9
        Optimal for 0.6m telescopes (Mag 10 to 17). Has true Johnson B and V.
        """
        mirrors = [
            'vizier.cfa.harvard.edu',  
            'vizier.cds.unistra.fr',   
            'vizier.nao.ac.jp',        
            'vizier.iucaa.in',         
            'vizier.china-vo.org'      
        ]
        
        for mirror in mirrors:
            print(f"  -> Attempting APASS connection via mirror: {mirror}")
            # FIX: Request all columns '*' to prevent strict naming crashes
            v = Vizier(columns=['RAJ2000', 'DEJ2000', '*'], row_limit=1000)
            v.TIMEOUT = 15  
            v.VIZIER_SERVER = mirror
            try:
                result = v.query_region(center_coord, radius=radius, catalog='II/336/apass9')
                if len(result) > 0:
                    return result[0]
            except Exception as e:
                print(f"  -> [WARNING] Mirror {mirror} blocked/failed: {e}")
                continue
        
        return None

    def _query_panstarrs(self, center_coord, radius):
        """
        PRIORITY 2 CATALOG: Pan-STARRS DR2
        Massive depth, but bright stars (>Mag 14) may be saturated. 
        Uses NASA MAST REST API to guarantee connection bypass of ISP blocks.
        """
        try:
            ra_deg = center_coord.ra.deg
            dec_deg = center_coord.dec.deg
            radius_deg = radius.to(u.deg).value

            url = f"https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv?ra={ra_deg}&dec={dec_deg}&radius={radius_deg}"
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            
            with urllib.request.urlopen(req, timeout=30) as response:
                csv_text = response.read().decode('utf-8')

            lines = [line.strip() for line in csv_text.strip().split('\n') if line.strip()]
            if len(lines) < 2: return None
                
            header = [h.strip() for h in lines[0].split(',')]
            ra_idx, dec_idx, r_idx = header.index('raMean'), header.index('decMean'), header.index('rMeanPSFMag')
            g_idx = header.index('gMeanPSFMag') if 'gMeanPSFMag' in header else -1
            i_idx = header.index('iMeanPSFMag') if 'iMeanPSFMag' in header else -1
            
            ra_list, dec_list, g_list, r_list, i_list = [], [], [], [], []
            
            for line in lines[1:]:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) <= max(ra_idx, dec_idx, r_idx): continue
                    
                try:
                    ra_list.append(float(parts[ra_idx]))
                    dec_list.append(float(parts[dec_idx]))
                except ValueError: continue
                
                def parse_mag(idx):
                    if idx == -1 or idx >= len(parts): return np.nan
                    try:
                        fval = float(parts[idx])
                        return fval if fval > -99.0 else np.nan
                    except ValueError: return np.nan
                        
                g_list.append(parse_mag(g_idx))
                r_list.append(parse_mag(r_idx))
                i_list.append(parse_mag(i_idx))
                
            mast_table = Table([ra_list, dec_list, g_list, r_list, i_list], names=('RAJ2000', 'DEJ2000', 'gmag', 'rmag', 'imag'))
            return mast_table
            
        except Exception:
            return None

    def _query_sdss(self, center_coord, radius):
        """
        PRIORITY 3 CATALOG: SDSS DR12
        Excellent fallback. Uses raw SQL to bypass Python radius limits.
        """
        try:
            ra_min = center_coord.ra.deg - radius.to(u.deg).value
            ra_max = center_coord.ra.deg + radius.to(u.deg).value
            dec_min = center_coord.dec.deg - radius.to(u.deg).value
            dec_max = center_coord.dec.deg + radius.to(u.deg).value
            
            sql = f"SELECT ra as RAJ2000, dec as DEJ2000, u, g, r, i, z FROM PhotoObj WHERE ra BETWEEN {ra_min} AND {ra_max} AND dec BETWEEN {dec_min} AND {dec_max} AND type=6"
            res = SDSS.query_sql(sql)
            return res
        except Exception:
            return None

    def _fetch_reference_catalog(self, center_coord, radius, filter_band):
        """
        The Priority Queue Manager.
        Decides which catalogs to query in what order based on the requested filter.
        """
        fb = filter_band.upper().replace('MAG', '').replace('-BAND', '').replace('_', '').strip()
        
        # Priority mapping: APASS is best for 0.6m scope, then Pan-STARRS, then SDSS.
        # FIX: Updated APASS names to perfectly match VizieR's apostrophe formatting.
        priority_queue = []
        if fb == 'R':
            priority_queue = [('APASS', "r'mag", self._query_apass), ('Pan-STARRS', 'rmag', self._query_panstarrs), ('SDSS', 'r', self._query_sdss)]
        elif fb == 'G':
            priority_queue = [('APASS', "g'mag", self._query_apass), ('Pan-STARRS', 'gmag', self._query_panstarrs), ('SDSS', 'g', self._query_sdss)]
        elif fb == 'I':
            priority_queue = [('APASS', "i'mag", self._query_apass), ('Pan-STARRS', 'imag', self._query_panstarrs), ('SDSS', 'i', self._query_sdss)]
        elif fb == 'V':
            priority_queue = [('APASS', 'Vmag', self._query_apass)] # Only APASS has true V
        elif fb == 'B':
            priority_queue = [('APASS', 'Bmag', self._query_apass)] # Only APASS has true B
        else:
            raise ValueError(f"Unsupported scientific filter band: {filter_band}")

        print(f"\n  -> Catalog Priority Sequence for Filter '{fb}':")
        for idx, (cat_name, col_name, _) in enumerate(priority_queue, 1):
            print(f"     {idx}. {cat_name} (Target Column: {col_name})")

        # Execute the fallback cascade
        for cat_name, col_name, query_func in priority_queue:
            print(f"  -> [ATTEMPTING] Fetching from {cat_name}...")
            
            raw_cat = query_func(center_coord, radius)
            
            if raw_cat is not None and len(raw_cat) > 0:
                # Find the correct magnitude column 
                actual_col = None
                if col_name in raw_cat.colnames:
                    actual_col = col_name
                else:
                    # FIX: Strip out apostrophes and underscores for bulletproof fuzzy matching
                    clean_target = col_name.replace("'", "").replace("_", "").lower()
                    possible_cols = [c for c in raw_cat.colnames if clean_target in c.replace("'", "").replace("_", "").lower() and 'mag' in c.lower()]
                    if possible_cols: actual_col = possible_cols[0]

                if actual_col:
                    # Clean out NaNs and masked values
                    valid_mask = ~np.isnan(raw_cat[actual_col])
                    if hasattr(raw_cat[actual_col], 'mask'):
                        valid_mask &= ~raw_cat[actual_col].mask
                        
                    clean_cat = raw_cat[valid_mask]
                    
                    if len(clean_cat) > 0:
                        print(f"  -> [SUCCESS] Extracted {len(clean_cat)} valid reference stars from {cat_name}.")
                        # Standardize the output so the math engine doesn't have to care which catalog won
                        std_table = Table()
                        std_table['RAJ2000'] = clean_cat['RAJ2000']
                        std_table['DEJ2000'] = clean_cat['DEJ2000']
                        std_table['Ref_Mag'] = clean_cat[actual_col]
                        return std_table
                    else:
                        print(f"  -> [FAILED] Stars found in {cat_name}, but none had valid {col_name} data.")
                else:
                    print(f"  -> [FAILED] Expected column '{col_name}' was missing from {cat_name} result. Found: {raw_cat.colnames}")
            else:
                print(f"  -> [FAILED] {cat_name} returned no stars or connection failed.")

        raise RuntimeError(f"All prioritized catalogs failed for filter '{fb}'. Cannot perform flux calibration.")

    def calculate_local_zero_point(self, master_science_fits, filter_band='rmag'):
        """
        The Swope Method: Finds all stars in the image, cross-matches them with the 
        winning catalog from the Priority Queue, and calculates a median Zero Point.
        """
        print("\n--- Calculating Local Zero Point (Differential Photometry) ---")
        with fits.open(master_science_fits) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header)
            
            ny, nx = data.shape
            center_coord = wcs.pixel_to_world(nx/2, ny/2)
            pixel_scale = np.abs(wcs.pixel_scale_matrix[0,0]) * u.deg
            search_radius = (nx/2) * pixel_scale
            
            # Fetch the best available catalog using the Fallback Queue
            catalog_table = self._fetch_reference_catalog(center_coord, search_radius, filter_band)
            catalog_coords = SkyCoord(ra=catalog_table['RAJ2000'], dec=catalog_table['DEJ2000'], unit=(u.deg, u.deg))
            
            mean, median, std = sigma_clipped_stats(data, sigma=3.0)
            finder = DAOStarFinder(fwhm=self.fwhm, threshold=self.threshold * std)
            detected_sources = finder.find_stars(data - median)
            
            detected_coords = wcs.pixel_to_world(detected_sources['xcentroid'], detected_sources['ycentroid'])
            idx, d2d, _ = match_coordinates_sky(detected_coords, catalog_coords)
            match_mask = d2d < 2.0 * u.arcsec
            
            calculated_zps = []
            
            for i in np.where(match_mask)[0]:
                cat_idx = idx[i]
                cat_mag = catalog_table['Ref_Mag'][cat_idx]
                
                position = (detected_sources['xcentroid'][i], detected_sources['ycentroid'][i])
                aperture = CircularAperture(position, r=self.fwhm * 2)
                annulus = CircularAnnulus(position, r_in=self.fwhm * 3, r_out=self.fwhm * 4)
                
                phot_table = aperture_photometry(data, [aperture, annulus])
                bkg_mean = phot_table['aperture_sum_1'][0] / annulus.area
                instrumental_flux = phot_table['aperture_sum_0'][0] - (bkg_mean * aperture.area)
                
                if instrumental_flux > 0:
                    inst_mag = -2.5 * np.log10(instrumental_flux)
                    star_zp = cat_mag - inst_mag
                    calculated_zps.append(star_zp)
            
            if len(calculated_zps) == 0:
                raise ValueError("Could not match any local stars with the catalog to determine a Zero Point.")
                
            self.zero_point = np.median(calculated_zps)
            print(f"  -> Successfully cross-matched {len(calculated_zps)} local stars.")
            print(f"  -> Local Differential Zero Point Locked: {self.zero_point:.4f} mag")
            
        return self.zero_point

    def export_flux_calibrated_image(self, target_fits, output_filename="calibrated_flux_science.fits"):
        if self.zero_point is None:
            raise RuntimeError("Zero Point missing. Run calculate_local_zero_point first.")
            
        print(f"\n--- Generating Flux-Calibrated FITS ({output_filename}) ---")
        conversion_factor = 3631.0 * (10 ** (-self.zero_point / 2.5))
        
        with fits.open(target_fits) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            calibrated_data = data * conversion_factor
            
            header['BUNIT'] = 'Jy'
            header['MAGZERO'] = self.zero_point
            header['FLUXCAL'] = conversion_factor
            header['HISTORY'] = f'Flux calibrated to Janskys using ZP: {self.zero_point:.4f}'
            
            fits.writeto(output_filename, calibrated_data, header, overwrite=True)
            print(f"  -> Success! Image pixels converted from electrons to Janskys.")

    def generate_full_catalog(self, target_fits, output_csv="object_catalog.csv", output_segmap="segmentation_map.fits"):
        if self.zero_point is None:
            raise RuntimeError("Zero Point missing. Run calculate_local_zero_point first.")
            
        print(f"\n--- Generating Full Object Catalog and Segmentation Map ---")
        with fits.open(target_fits) as hdul:
            data = hdul[0].data
            wcs = WCS(hdul[0].header)
            
            bkg_estimator = MedianBackground()
            bkg = Background2D(data, (50, 50), filter_size=(3, 3), bkg_estimator=bkg_estimator)
            data_sub = data - bkg.background
            threshold = self.threshold * bkg.background_rms
            
            kernel = Gaussian2DKernel(x_stddev=self.fwhm / 2.35)
            kernel.normalize()
            convolved_data = convolve(data_sub, kernel)
            
            segment_map = detect_sources(convolved_data, threshold, npixels=5)
            print(f"  -> Deblending touching objects...")
            segment_map = deblend_sources(data_sub, segment_map, npixels=5, nlevels=32, contrast=0.001)
            print(f"  -> After deblending: {segment_map.nlabels} unique objects (Stars + Galaxies).")
            
            if output_segmap:
                seg_header = hdul[0].header.copy()
                seg_header['BUNIT'] = 'ID'
                seg_header['HISTORY'] = 'Segmentation map generated by UniversalPhotometryEngine'
                fits.writeto(output_segmap, segment_map.data.astype(np.int32), seg_header, overwrite=True)
            
            cat = SourceCatalog(data_sub, segment_map, wcs=wcs)
            tbl = cat.to_table()
            
            valid_mask = tbl['segment_flux'] > 0
            tbl = tbl[valid_mask]
            cat = cat[valid_mask] 
            
            final_catalog = Table()
            final_catalog['ID'] = tbl['label']
            final_catalog['RA_deg'] = cat.sky_centroid.ra.deg
            final_catalog['Dec_deg'] = cat.sky_centroid.dec.deg
            final_catalog['X_pix'] = tbl['xcentroid']
            final_catalog['Y_pix'] = tbl['ycentroid']
            final_catalog['Instrumental_Flux'] = tbl['segment_flux']
            final_catalog['Absolute_Mag'] = -2.5 * np.log10(tbl['segment_flux']) + self.zero_point
            final_catalog['Area_pixels'] = tbl['area']
            
            if 'ellipticity' in tbl.colnames: final_catalog['Ellipticity'] = tbl['ellipticity']
            elif hasattr(cat, 'ellipticity'): final_catalog['Ellipticity'] = cat.ellipticity
            elif 'eccentricity' in tbl.colnames: final_catalog['Ellipticity'] = tbl['eccentricity']
            elif hasattr(cat, 'eccentricity'): final_catalog['Ellipticity'] = cat.eccentricity
            else: final_catalog['Ellipticity'] = np.nan
            
            final_catalog.write(output_csv, format='csv', overwrite=True)
            print(f"  -> Catalog saved successfully to {output_csv}!")