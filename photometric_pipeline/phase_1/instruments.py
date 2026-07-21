# instruments.py
from instrument_base import InstrumentProfile

class ITelescopeNetworkProfile(InstrumentProfile):
    def __init__(self):
        # Fallback database in case headers are missing values.
        self.hardware_database = {
            'T32': {'gain': 1.0, 'read_noise': 9.0, 'fringes': False},
            'T24': {'gain': 1.4, 'read_noise': 7.0, 'fringes': False},
            'T68': {'gain': 0.77, 'read_noise': 1.5, 'fringes': False}, 
            'T11': {'gain': 1.3, 'read_noise': 8.0, 'fringes': True},   
            'DEFAULT': {'gain': 1.0, 'read_noise': 10.0, 'fringes': False}
        }

    @property
    def name(self):
        return "iTelescope Universal Network"

    def _get_telescope_id(self, header):
        """Helper function to extract which iTelescope took the image."""
        tel_string = header.get('TELESCOP', header.get('INSTRUME', 'UNKNOWN')).upper()
        
        for tel_id in self.hardware_database.keys():
            if tel_id in tel_string:
                return tel_id
        return 'DEFAULT'

    def get_amplifiers(self, hdul):
        return [hdul[0]]

    def get_image_type(self, header):
        img_type = header.get('IMAGETYP', '').lower()
        if 'bias' in img_type: return 'bias'
        if 'dark' in img_type: return 'dark'
        if 'flat' in img_type: return 'flat'
        return 'science'

    def get_exposure(self, header):
        return header.get('EXPTIME', 0.0)

    def get_gain(self, header):
        # Safely attempt to extract gain
        for key in ['EGAIN', 'GAIN', 'SYSGAIN']:
            val = header.get(key)
            if val is not None:
                try:
                    return float(val)
                except ValueError:
                    continue # If it's weird text, ignore and keep looking
                    
        tel_id = self._get_telescope_id(header)
        return self.hardware_database[tel_id]['gain']

    def get_read_noise(self, header):
        # Safely attempt to extract read noise
        for key in ['READNOIS', 'RDNOISE', 'E-NOISE']:
            val = header.get(key)
            if val is not None:
                try:
                    return float(val)
                except ValueError:
                    continue # If it's a string like 'Mode0', ignore and keep looking
                    
        tel_id = self._get_telescope_id(header)
        return self.hardware_database[tel_id]['read_noise']

    def get_overscan_region(self, header):
        return None 
        
    def needs_fringe_correction(self, header):
        tel_id = self._get_telescope_id(header)
        telescope_fringes = self.hardware_database[tel_id]['fringes']
        
        filter_name = header.get('FILTER', '').lower()
        fringe_filters = ['z', 'y', 'ha', 'sii', 'luminance'] 
        
        if telescope_fringes and any(f in filter_name for f in fringe_filters):
            return True
            
        return False

    def get_filter(self, header):
        """Extracts the filter name and standardizes it."""
        return header.get('FILTER', 'UNKNOWN').strip().upper()