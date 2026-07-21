# instrument_base.py

class InstrumentProfile:
    """Base class for mapping telescope-specific FITS data to a standard model."""
    
    @property
    def name(self):
        return "Generic Instrument"

    def get_amplifiers(self, hdul):
        """Returns a list of FITS extensions that contain actual science/calibration data."""
        raise NotImplementedError

    def get_image_type(self, header):
        """Standardizes image types to: 'bias', 'dark', 'flat', or 'science'."""
        raise NotImplementedError

    def get_exposure(self, header):
        raise NotImplementedError

    def get_gain(self, header):
        raise NotImplementedError

    def get_read_noise(self, header):
        raise NotImplementedError
    
    def get_filter(self, header):
        """Extracts the filter name and standardizes it to uppercase."""
        raise NotImplementedError

    def get_overscan_region(self, header):
        """Returns the overscan slice (e.g., '[1:2048, 2049:2100]') or None."""
        return None
        
    def needs_fringe_correction(self, header):
        return False