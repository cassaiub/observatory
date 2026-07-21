# run_phase2.py
import os
import sys
import glob
import logging
import argparse
import gc
import numpy as np
import numpy.ma as ma
from astropy.stats import sigma_clip, mad_std
from scipy.ndimage import median_filter
from skimage.transform import warp
import astroalign as aa

from phase2_models import TargetGroup
from phase2_hardware import HardwareManager
from phase2_io import FITSHandler
from phase2_math import MathEngine
from phase2_wcs import WCSSolver
from phase2_visuals import VisualQAGenerator

import warnings
warnings.filterwarnings('ignore')

class IntegrationPipeline:
    def __init__(self, input_dir, keep_temps=False):
        self.input_dir = input_dir
        self.keep_temps = keep_temps
        self.hw = HardwareManager()
        self.target_groups = []
        self.run_dir = None
        self.logger = None

    def setup(self):
        cores = self.hw.allocate_resources()
        
        files = [f for f in glob.glob(os.path.join(self.input_dir, "*.fit*")) if os.path.dirname(f) == self.input_dir]
        buckets = {}
        t_set, f_set, s_set = set(), set(), set()
        
        print("\n[~] Skimming headers and file sizes...")
        total_size_mb = 0
        valid_files = 0
        
        for f in files:
            if "Master_" in os.path.basename(f): continue
            meta = FITSHandler.extract_metadata(f)
            if not meta: continue
            
            key = f"{meta['object']}_{meta['filter']}_{meta['hardware_id']}_{meta['exposure']}s"
            if key not in buckets: 
                buckets[key] = TargetGroup(key, meta)
            buckets[key].raw_files.append(f)
            t_set.add(meta['object']); f_set.add(meta['filter']); s_set.add(meta['hardware_id'].split('_')[0])
            
            total_size_mb += os.path.getsize(f) / (1024 * 1024)
            valid_files += 1
            
        self.target_groups = list(buckets.values())
        if not self.target_groups:
            print("[ERROR] No valid FITS found."); sys.exit(1)

        # --- RESTORED: Resource Estimation ---
        max_stack = max(len(g.raw_files) for g in self.target_groups) if self.target_groups else 0
        self.hw.print_estimations(valid_files, max_stack, total_size_mb)

        base_name = f"{'-'.join(sorted(t_set))}_{'-'.join(sorted(f_set))}_{'-'.join(sorted(s_set))}"
        count = 1
        while True:
            self.run_dir = os.path.join(self.input_dir, f"Run{count:02d}_{base_name}")
            if not os.path.exists(self.run_dir): os.makedirs(self.run_dir); break
            count += 1

        self.logger = logging.getLogger('Phase2_OO')
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logging.FileHandler(os.path.join(self.run_dir, "phase2_pipeline.log")))
        self.logger.addHandler(logging.StreamHandler())
        
        self.logger.info("==================================================")
        self.logger.info(f"PHASE 2: UNIFIED PIPELINE (Stacks + WCS + QA) - Cores: {cores}")
        self.logger.info("==================================================")

    def execute(self):
        wcs_solver = WCSSolver(self.logger)
        failed_wcs = []

        for group in self.target_groups:
            tot = len(group.raw_files)
            if tot < 2: continue
                
            self.logger.info(f"\n--- Stacking: {group.group_key} ({tot} frames) ---")
            
            best_snr, best_noise, group.anchor_filepath = 0, float('inf'), group.raw_files[0]
            for f in group.raw_files:
                td = FITSHandler.load_data(f)
                tc, _ = MathEngine.extract_2d_background(td)
                n = mad_std(tc, ignore_nan=True)
                snr = np.nanpercentile(tc, 99) / n if n > 0 else 0
                if snr > best_snr: best_snr, best_noise, group.anchor_filepath = snr, n, f
            
            self.logger.info(f"[+] Anchor set: {os.path.basename(group.anchor_filepath)}")
            
            r_data, _ = MathEngine.extract_2d_background(FITSHandler.load_data(group.anchor_filepath))
            stack = [r_data]
            weights = [1.0 / (best_noise**2) if best_noise > 0 else 1.0]
            
            # --- RESTORED: QA Tracking Arrays ---
            qa_stars, qa_scales, qa_noises = [], [], [best_noise]
            
            for f in group.raw_files:
                if f == group.anchor_filepath: continue
                t_data, _ = MathEngine.extract_2d_background(FITSHandler.load_data(f))
                
                try:
                    tf, (src, dst) = aa.find_transform(median_filter(t_data, 3), median_filter(r_data, 3), detection_sigma=1.5, min_area=4)
                    scale = MathEngine.calc_scale(t_data, r_data, src, dst)
                    if not (0.65 <= scale <= 1.5):
                        self.logger.warning(f"  [!] REJECTED: Scale {scale:.2f}x")
                        continue
                        
                    reg_data = warp(t_data, inverse_map=tf.inverse, output_shape=r_data.shape, order=3, cval=np.nan, preserve_range=True)
                    norm_data = reg_data * scale
                    fn = mad_std(norm_data, ignore_nan=True)
                    
                    stack.append(norm_data)
                    weights.append(1.0 / (fn**2) if fn > 0 else 0)
                    
                    # Track QA
                    qa_stars.append(len(src))
                    qa_scales.append(scale)
                    qa_noises.append(fn)
                    
                    self.logger.info(f"  [+] Aligned: {os.path.basename(f)} | Stars: {len(src)} | Scale: {scale:.2f}x")
                except Exception as e:
                    self.logger.error(f"  [-] Alignment failed for {os.path.basename(f)}: {e}")

            group.successful_frames = len(stack)
            if group.successful_frames < 3: continue

            self.logger.info("  [*] Integrating cube...")
            cube, w_arr = np.array(stack), np.array(weights) / sum(weights)
            
            if group.successful_frames < 5:
                mask = (cube == np.nanmin(cube, axis=0)) | (cube == np.nanmax(cube, axis=0)) | np.isnan(cube)
                master = ma.average(ma.masked_array(cube, mask=mask), axis=0, weights=w_arr).filled(np.nan)
            else:
                master = ma.average(sigma_clip(cube, sigma=3.0, maxiters=3, axis=0, masked=True), axis=0, weights=w_arr).filled(np.nan)
                
            master += 100.0
            
            group.master_filepath = os.path.join(self.run_dir, f"Master_{group.group_key}_{group.successful_frames}fr.fits")
            FITSHandler.save_master(master, group.anchor_filepath, group.master_filepath, group.successful_frames, group.meta)
            
            # --- RESTORED: Automated QA Telemetry ---
            self.logger.info(f"\n   >>> AUTOMATED QA TELEMETRY <<<")
            if qa_stars: self.logger.info(f"   [*] STAR MATCH : Avg: {int(np.mean(qa_stars))} | Min: {int(np.min(qa_stars))}")
            if qa_scales: self.logger.info(f"   [*] FLUX SCALE : Avg: {np.mean(qa_scales):.3f}x")

            c_anchor = MathEngine.center_crop(r_data)
            c_master = MathEngine.center_crop(master)
            
            anchor_noise = mad_std(c_anchor, ignore_nan=True)
            master_noise = mad_std(c_master, ignore_nan=True)
            
            actual_boost = anchor_noise / master_noise if master_noise > 0 else 0
            variances = np.array(qa_noises) ** 2
            expected_noise = np.sqrt(1.0 / np.sum(1.0 / variances))
            weighted_max = anchor_noise / expected_noise if expected_noise > 0 else 0
            efficiency = (actual_boost / weighted_max) * 100 if weighted_max > 0 else 0
            
            self.logger.info(f"   [*] PERFORMANCE: SNR Boost: {actual_boost:.2f}x | Efficiency: {efficiency:.1f}%")
            if efficiency < 50.0:
                self.logger.warning("       -> WARNING: Integration efficiency critically low (<50%).")
            self.logger.info("   --------------------------------------\n")

            del stack, cube, master, c_anchor, c_master; gc.collect()

            self.logger.info("  [*] Initiating Astrometry WCS Solve...")
            if not wcs_solver.solve(group, keep_temps=self.keep_temps):
                failed_wcs.append(group)

        if failed_wcs and (wcs_solver.global_anchor_ra is not None):
            self.logger.info("\n--- INITIATING WCS RESCUE PASS ---")
            for g in failed_wcs:
                wcs_solver.solve(g, is_rescue=True, keep_temps=self.keep_temps)

        VisualQAGenerator.generate_pdf(self.target_groups, self.run_dir, self.logger)
        self.logger.info("\n==================================================")
        self.logger.info("PHASE 2 PIPELINE COMPLETE.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified Phase 2 OO Pipeline")
    parser.add_argument("data_dir", type=str, help="Directory containing calibrated FITS")
    parser.add_argument("--keep-temps", action="store_true", help="Keep temporary files generated by solve-field.") # RESTORED
    args = parser.parse_args()

    if not os.path.exists(args.data_dir):
        print(f"Directory {args.data_dir} does not exist.")
        sys.exit(1)

    pipeline = IntegrationPipeline(os.path.abspath(args.data_dir), keep_temps=args.keep_temps)
    pipeline.setup()
    
    proceed = input("\nProceed with Unified Integration? [Y/n]: ").strip().lower()
    if proceed in ['y', '']: pipeline.execute()