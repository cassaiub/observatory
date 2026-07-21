# data_models.py
from astropy.io import fits
from astropy.nddata import CCDData
import astropy.units as u

class StandardCCD:
    """A standardized wrapper containing the AstroPy CCDData and unified metadata."""
    def __init__(self, ccd_data, meta):
        self.ccd = ccd_data
        self.meta = meta

def load_standardized_ccds(filepath, instrument):
    """
    Reads a FITS file using the provided instrument profile.
    ALWAYS returns a list of StandardCCD objects (length 1 for amateur, length N for MEF).
    """
    standardized_list = []
    
    with fits.open(filepath) as hdul:
        global_header = hdul[0].header
        amplifiers = instrument.get_amplifiers(hdul)
        
        for amp in amplifiers:
            # Merge global header with extension header for complete metadata
            header = global_header.copy()
            header.update(amp.header)
            
            # Extract standard metadata using the instrument driver
            meta = {
                'image_type': instrument.get_image_type(header),
                'exposure': instrument.get_exposure(header),
                'gain': instrument.get_gain(header),
                'read_noise': instrument.get_read_noise(header),
                'overscan': instrument.get_overscan_region(header),
                'fringe_needed': instrument.needs_fringe_correction(header),
                'filter': instrument.get_filter(header) # NEW ADDITION
            }
            
            # Create the AstroPy CCD object
            ccd = CCDData(amp.data, meta=header, unit=u.adu)
            standardized_list.append(StandardCCD(ccd, meta))
            
    return standardized_list