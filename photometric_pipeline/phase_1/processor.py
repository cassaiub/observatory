# processor.py
import ccdproc
import numpy as np
import astropy.units as u
import astroscrappy

class UniversalProcessor:
    def __init__(self, master_bias=None, master_dark=None, master_flats=None, bpm=None):
        self.master_bias = master_bias
        self.master_dark = master_dark
        self.master_flats = master_flats or {} 
        self.bpm = bpm

    def process_science_frame(self, standard_ccd_list):
        processed_ccds = []
        
        for std_ccd in standard_ccd_list:
            ccd = std_ccd.ccd
            meta = std_ccd.meta
            
            # 1. Overscan Subtraction
            if meta['overscan']:
                ccd = ccdproc.subtract_overscan(ccd, overscan=ccd[:, meta['overscan']], median=True)
                ccd = ccdproc.trim_image(ccd, fits_section=meta['overscan'])
                
            # 2. Bias Subtraction
            if self.master_bias:
                # Safety check for binning mismatch
                if ccd.shape != self.master_bias.shape:
                    print(f"  -> [ERROR] Dimension mismatch! Science frame is {ccd.shape}, but Master Bias is {self.master_bias.shape}.")
                    print(f"  -> Check your camera Binning settings. Skipping this frame.")
                    continue
                ccd = ccdproc.subtract_bias(ccd, self.master_bias)
                
            # 3. Dark Subtraction
            if self.master_dark:
                # We feed it the exact numeric value we extracted, converted to an Astropy Time unit.
                # This bypasses all FITS header lookups and satisfies the ccdproc library requirements.
                ccd = ccdproc.subtract_dark(
                    ccd, 
                    self.master_dark, 
                    dark_exposure=meta['exposure'] * u.s, 
                    data_exposure=meta['exposure'] * u.s, 
                    scale=False
                )
                
            # 4. Flat Fielding
            filter_name = meta['filter']
            matching_flat = self.master_flats.get(filter_name)
            
            if matching_flat:
                ccd = ccdproc.flat_correct(ccd, matching_flat)
            else:
                print(f"  -> [Warning] No Master Flat found for filter: {filter_name}. Skipping flat fielding.")
                
            # 5. Bad Pixel Masking
            if self.bpm is not None:
                ccd.mask = self.bpm
                
            # 6. Cosmic Ray Rejection
            crmask, clean_data = astroscrappy.detect_cosmics(
                ccd.data, inmask=ccd.mask, gain=meta['gain'], readnoise=meta['read_noise'],
                sigclip=4.5, sigfrac=0.3, objlim=5.0
            )
            ccd.data = clean_data
            if ccd.mask is not None:
                ccd.mask = np.logical_or(ccd.mask, crmask)
            else:
                ccd.mask = crmask
                
            # 7. Gain Correction (ADU -> Electrons)
            ccd = ccdproc.gain_correct(ccd, meta['gain'] * (u.electron / u.adu))
            
            if meta['fringe_needed']:
                pass 
                
            processed_ccds.append(ccd)
            
        return processed_ccds