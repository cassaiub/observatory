# phase2_io.py
import numpy as np
from astropy.io import fits

class FITSHandler:
    """Handles FITS reading, taxonomy mapping, and Master saving."""
    FILTER_MAP = {
        'RED': 'Red', 'GREEN': 'Green', 'BLUE': 'Blue',
        'LUM': 'L', 'LUMINANCE': 'L', 'CLEAR': 'L',
        'R': 'R_Photo', 'V': 'V_Photo', 'B': 'B_Photo', 'U': 'U_Photo', 'I': 'I_Photo',
        'R-PHOTOMETRIC': 'R_Photo', 'V-PHOTOMETRIC': 'V_Photo',
        'HALPHA': 'Ha', 'H-ALPHA': 'Ha', 'HA': 'Ha',
        'HA-7NM': 'Ha_7nm', 'HA-3NM': 'Ha_3nm',
        'S-II': 'SII', 'SII': 'SII', 'S2': 'SII',
        'O-III': 'OIII', 'OIII': 'OIII', 'O3': 'OIII'
    }

    @staticmethod
    def extract_metadata(filepath):
        try:
            with fits.open(filepath) as hdul:
                header = hdul[0].header
                meta = {}
                meta['object'] = header.get('OBJECT', header.get('TARGET', 'Unknown')).replace(" ", "").upper()
                
                raw_filter = header.get('FILTER', 'NoFilter').strip().upper()
                bandwidth = header.get('BANDWID', '')
                hw_suffix = f"_{bandwidth}nm" if bandwidth else ""
                
                if "NM" in raw_filter and raw_filter not in FITSHandler.FILTER_MAP:
                    meta['filter'] = raw_filter.replace("-", "_")
                else:
                    meta['filter'] = FITSHandler.FILTER_MAP.get(raw_filter, raw_filter) + hw_suffix
                
                meta['hardware_id'] = header.get('INSTRUME', 'Unknown_Cam').replace(" ", "")
                raw_exp = float(header.get('EXPTIME', header.get('EXPOSURE', 0.0)))
                meta['exposure'] = int(5 * round(raw_exp / 5))
                meta['pixel_size_um'] = float(header.get('XPIXSZ', 0.0))
                meta['focal_length_mm'] = float(header.get('FOCALLEN', 0.0))
                meta['binning'] = int(header.get('XBINNING', 1))
                return meta
        except Exception:
            return None

    @staticmethod
    def load_data(filepath):
        with fits.open(filepath) as hdul:
            return hdul[0].data.astype(np.float32)

    @staticmethod
    def save_master(master_data, ref_filepath, output_filename, num_frames, meta):
        with fits.open(ref_filepath) as hdul:
            master_header = hdul[0].header
            
        pixel_scale = 0.0
        if meta['focal_length_mm'] > 0 and meta['pixel_size_um'] > 0:
            pixel_scale = (meta['pixel_size_um'] / meta['focal_length_mm']) * 206.265
            
        master_header['STACKCNT'] = (num_frames, 'Total frames integrated')
        master_header['TOT_EXP']  = (num_frames * meta['exposure'], '[s] Total integration time')
        master_header['PIXSCALE'] = (round(pixel_scale, 4), '[arcsec/pix] Plate scale')
        
        bg_val = float(np.nanmedian(master_data))
        master_header['BKGND'] = (bg_val if not np.isnan(bg_val) else 0.0, 'Median background')
        master_header['HISTORY'] = 'PHASE 2: SNR Anchor, Rejection, Weighted, 2D BG Subtracted'
        
        clean_data = np.nan_to_num(master_data, nan=0.0)
        hdu = fits.PrimaryHDU(clean_data, header=master_header)
        fits.HDUList([hdu]).writeto(output_filename, overwrite=True)
        return pixel_scale