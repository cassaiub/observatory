# phase2_wcs.py
import os
import subprocess
import shutil
from astropy.io import fits

class WCSSolver:
    def __init__(self, logger):
        self.logger = logger
        
        # 1. Find the absolute path of the directory containing THIS script
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        
        # 2. Build paths relative to this script's location
        # NOTE: If your index files on the HPC are elsewhere, you must move them here 
        # or change this back to the absolute HPC path!
        self.index_dir = os.path.join(self.base_dir, "astrometry_data")
        self.config_path = os.path.join(self.base_dir, "astrometry.cfg")
        
        self.global_anchor_ra = None
        self.global_anchor_dec = None
        
        self._ensure_config()

    def _ensure_config(self):
        """Creates or patches the astrometry.cfg file dynamically."""
        if not os.path.exists(self.config_path):
            self.logger.info(f"[*] Generating local Astrometry config at {self.config_path}")
            try:
                with open(self.config_path, 'w') as f:
                    f.write(f"add_path {self.index_dir}\n")
                    f.write("autoindex\n")
                    f.write("inparallel\n")
            except Exception as e:
                self.logger.error(f"[!] Failed to create config file: {e}")
            return

        try:
            with open(self.config_path, 'r') as f: 
                content = f.read()
            if self.index_dir not in content:
                self.logger.info("[*] Patching existing astrometry.cfg with local index...")
                with open(self.config_path, 'a') as f:
                    f.write(f"\n# Phase 2 Pipeline Auto-Patch\nadd_path {self.index_dir}\nautoindex\ninparallel\n")
        except Exception as e:
            self.logger.error(f"[!] Failed to patch config file: {e}")

    def solve(self, group, is_rescue=False, keep_temps=False):
        filepath = group.master_filepath
        if not filepath or not os.path.exists(filepath): return False

        output_fits = os.path.splitext(filepath)[0] + "_wcs.fits"
        flag_file = os.path.splitext(filepath)[0] + ".solved"
        for f in [output_fits, flag_file]:
            if os.path.exists(f): os.remove(f)

        scale_args = []
        pix_scale = 0.0
        try:
            with fits.open(filepath) as hdul:
                hdr = hdul[0].header
                pix_scale = float(hdr.get('PIXSCALE', 0.0))
                ra_hint = self.global_anchor_ra or hdr.get('OBJCTRA') or hdr.get('RA')
                dec_hint = self.global_anchor_dec or hdr.get('OBJCTDEC') or hdr.get('DEC')

            if pix_scale > 0:
                scale_args.extend(["--scale-units", "arcsecperpix", "--scale-low", str(pix_scale*0.75), "--scale-high", str(pix_scale*1.25)])
            if ra_hint and dec_hint:
                ra_f, dec_f = str(ra_hint).strip().replace(" ", ":"), str(dec_hint).strip().replace(" ", ":")
                scale_args.extend(["--ra", ra_f, "--dec", dec_f, "--radius", "3.0"])
        except Exception: pass

        if pix_scale > 1.0: ds, obj, sig = "4", "150", "10"
        elif 0.0 < pix_scale <= 1.0:
            if any(nb in filepath for nb in ['_Ha_', '_OIII_', '_SII_']): ds, obj, sig = "1", "1000", "3"
            else: ds, obj, sig = "2", "600", "8"
        else: ds, obj, sig = "2", "300", "8" 

        self.logger.info(f"    -> WCS Profile: DS {ds}x | Objs: {obj} | Sigma: {sig}")
        
        cmd = ["solve-field", filepath, "--config", self.config_path, "--overwrite", "--no-plots", "--downsample", ds, "--objs", obj, "--sigma", sig, "--cpulimit", "60", "--new-fits", output_fits] + scale_args
        
        try:
            # Added text=True to capture the raw strings from the console
            res_fast = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            if not (os.path.exists(flag_file) and os.path.exists(output_fits)) and not is_rescue:
                self.logger.info("    [!] Fast solve failed. Engaging Blind Solve (300s)...")
                
                # LOG THE FAST SOLVE FAILURE
                self.logger.error("        --- ASTROMETRY FAST SOLVE LOG ---")
                if res_fast.stdout: self.logger.error(res_fast.stdout.strip())
                if res_fast.stderr: self.logger.error(res_fast.stderr.strip())
                self.logger.error("        ---------------------------------")
                
                cmd_blind = ["solve-field", filepath, "--config", self.config_path, "--overwrite", "--no-plots", "--downsample", ds, "--objs", obj, "--sigma", sig, "--cpulimit", "300", "--new-fits", output_fits]
                res_blind = subprocess.run(cmd_blind, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                # LOG THE BLIND SOLVE FAILURE
                if not (os.path.exists(flag_file) and os.path.exists(output_fits)):
                    self.logger.error("        --- ASTROMETRY BLIND SOLVE LOG ---")
                    if res_blind.stdout: self.logger.error(res_blind.stdout.strip())
                    if res_blind.stderr: self.logger.error(res_blind.stderr.strip())
                    self.logger.error("        ----------------------------------")
                    
        except FileNotFoundError:
            self.logger.error("    [-] FATAL WCS ERROR: 'solve-field' command not found.")
            self.logger.error("    [-] Ensure Astrometry.net is installed and in your system PATH.")
            return False

        if os.path.exists(flag_file) and os.path.exists(output_fits):
            with fits.open(output_fits, mode='update') as hdul:
                self.global_anchor_ra = hdul[0].header.get('CRVAL1')
                self.global_anchor_dec = hdul[0].header.get('CRVAL2')
                hdul[0].header['HISTORY'] = 'PHASE 2: WCS solved via Astrometry.net'
            shutil.move(output_fits, filepath)
            
            if not keep_temps: 
                self._clean_temps(filepath)
                
            group.wcs_solved = True
            self.logger.info("    [+] SUCCESS: WCS mapped.")
            return True
        else:
            self.logger.error("    [-] WCS FAILED.")
            if not keep_temps: 
                self._clean_temps(filepath)
            return False

    def _clean_temps(self, base_filepath):
        base = os.path.splitext(base_filepath)[0]
        for ext in ['-indx.xyls', '.axy', '.corr', '.match', '.rdls', '.solved', '.wcs']:
            try: os.remove(f"{base}{ext}")
            except OSError: pass
