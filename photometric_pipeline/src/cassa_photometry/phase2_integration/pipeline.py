"""Phase 2 orchestrator: align, stack, and astrometrically solve calibrated frames.

The stack is inverse-variance weighted and now carries a propagated ``ERR``
plane through to the master (previously the variance was computed only for a QA
log line and discarded). The arbitrary ``+100`` pedestal has been removed so the
output remains on a physical electron scale for photometry.
"""

import os
import gc
import glob
import warnings

import numpy as np
from astropy.stats import mad_std
from scipy.ndimage import median_filter
import astroalign as aa

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.fits_utils import build_dq
from cassa_photometry.phase2_integration.models import TargetGroup
from cassa_photometry.phase2_integration.hardware import HardwareManager
from cassa_photometry.phase2_integration.io import FITSHandler
from cassa_photometry.phase2_integration.math_utils import MathEngine
from cassa_photometry.phase2_integration.wcs import WCSSolver
from cassa_photometry.phase2_integration.visuals import VisualQAGenerator

warnings.filterwarnings("ignore")


def _variance_from(err, shape, fallback_noise):
    """Return a per-pixel variance plane, falling back to a flat scalar variance."""
    if err is not None:
        return np.asarray(err, dtype=np.float64) ** 2
    return np.full(shape, float(fallback_noise) ** 2, dtype=np.float64)


class IntegrationPipeline:
    def __init__(self, input_dir, keep_temps=False, config=None, logger=None):
        self.input_dir = input_dir
        self.keep_temps = keep_temps
        self.config = config or load_config()
        self.hw = HardwareManager()
        self.target_groups = []
        self.run_dir = None
        self.logger = logger

    def setup(self):
        cores = self.hw.allocate_resources()

        files = [f for f in glob.glob(os.path.join(self.input_dir, "*.fit*"))
                 if os.path.dirname(f) == self.input_dir]
        buckets = {}
        t_set, f_set, s_set = set(), set(), set()

        print("\n[~] Skimming headers and file sizes...")
        total_size_mb = 0
        valid_files = 0

        for f in files:
            if "Master_" in os.path.basename(f):
                continue
            meta = FITSHandler.extract_metadata(f)
            if not meta:
                continue

            key = f"{meta['object']}_{meta['filter']}_{meta['hardware_id']}_{meta['exposure']}s"
            buckets.setdefault(key, TargetGroup(key, meta)).raw_files.append(f)
            t_set.add(meta["object"])
            f_set.add(meta["filter"])
            s_set.add(meta["hardware_id"].split("_")[0])

            total_size_mb += os.path.getsize(f) / (1024 * 1024)
            valid_files += 1

        self.target_groups = list(buckets.values())
        if not self.target_groups:
            raise SystemExit("[ERROR] No valid FITS found.")

        max_stack = max(len(g.raw_files) for g in self.target_groups)
        self.hw.print_estimations(valid_files, max_stack, total_size_mb)

        base_name = f"{'-'.join(sorted(t_set))}_{'-'.join(sorted(f_set))}_{'-'.join(sorted(s_set))}"
        count = 1
        while True:
            self.run_dir = os.path.join(self.input_dir, f"Run{count:02d}_{base_name}")
            if not os.path.exists(self.run_dir):
                os.makedirs(self.run_dir)
                break
            count += 1

        self.logger = get_logger("cassa_integrate", run_dir=self.run_dir)
        self.logger.info("=" * 50)
        self.logger.info(f"PHASE 2: UNIFIED PIPELINE (Stacks + WCS + QA) - Cores: {cores}")
        self.logger.info("=" * 50)

    def execute(self):
        cfg = self.config.phase2
        wcs_solver = WCSSolver(self.logger, self.config)
        failed_wcs = []

        for group in self.target_groups:
            tot = len(group.raw_files)
            if tot < 2:
                continue

            self.logger.info(f"\n--- Stacking: {group.group_key} ({tot} frames) ---")

            # Anchor = highest-SNR frame.
            best_snr, best_noise, group.anchor_filepath = 0, float("inf"), group.raw_files[0]
            for f in group.raw_files:
                td = FITSHandler.load_data(f)
                tc, _ = MathEngine.extract_2d_background(td, cfg.background_box, cfg.background_filter)
                n = mad_std(tc, ignore_nan=True)
                snr = np.nanpercentile(tc, 99) / n if n > 0 else 0
                if snr > best_snr:
                    best_snr, best_noise, group.anchor_filepath = snr, n, f

            self.logger.info(f"[+] Anchor set: {os.path.basename(group.anchor_filepath)}")

            a_sci, a_err, _ = FITSHandler.load_planes(group.anchor_filepath)
            r_data, _ = MathEngine.extract_2d_background(a_sci, cfg.background_box, cfg.background_filter)
            r_var = _variance_from(a_err, r_data.shape, best_noise)

            stack = [r_data]
            var_stack = [r_var]
            weights = [1.0 / (best_noise ** 2) if best_noise > 0 else 1.0]
            qa_stars, qa_scales, qa_noises = [], [], [best_noise]

            for f in group.raw_files:
                if f == group.anchor_filepath:
                    continue
                t_sci, t_err, _ = FITSHandler.load_planes(f)
                t_data, _ = MathEngine.extract_2d_background(t_sci, cfg.background_box, cfg.background_filter)

                try:
                    tf, (src, dst) = aa.find_transform(
                        median_filter(t_data, 3), median_filter(r_data, 3),
                        detection_sigma=cfg.align_detection_sigma, min_area=cfg.align_min_area,
                    )
                    scale = MathEngine.calc_scale(t_data, r_data, src, dst)
                    if not (cfg.scale_min <= scale <= cfg.scale_max):
                        self.logger.warning(f"  [!] REJECTED: Scale {scale:.2f}x")
                        continue

                    reg_data = MathEngine.register_plane(t_data, tf.inverse, r_data.shape, cfg.warp_order)
                    norm_data = reg_data * scale
                    fn = mad_std(norm_data, ignore_nan=True)

                    t_var = _variance_from(t_err, t_data.shape, fn)
                    reg_var = MathEngine.register_variance(t_var, tf.inverse, r_data.shape)
                    norm_var = reg_var * (scale ** 2)

                    stack.append(norm_data)
                    var_stack.append(norm_var)
                    weights.append(1.0 / (fn ** 2) if fn > 0 else 0)

                    qa_stars.append(len(src))
                    qa_scales.append(scale)
                    qa_noises.append(fn)
                    self.logger.info(f"  [+] Aligned: {os.path.basename(f)} | Stars: {len(src)} | Scale: {scale:.2f}x")
                except Exception as exc:
                    self.logger.error(f"  [-] Alignment failed for {os.path.basename(f)}: {exc}")

            group.successful_frames = len(stack)
            if group.successful_frames < 3:
                continue

            self.logger.info("  [*] Integrating cube (inverse-variance weighted)...")
            cube = np.array(stack)
            varcube = np.array(var_stack)
            w_arr = np.array(weights) / sum(weights)
            use_clip = group.successful_frames >= cfg.sigma_clip_min_frames

            master, master_var = MathEngine.weighted_stack(
                cube, varcube, w_arr,
                sigma=cfg.stack_sigma, maxiters=cfg.stack_maxiters, use_sigma_clip=use_clip,
            )
            master_err = np.sqrt(master_var)
            master_dq = build_dq(master.shape, no_data=~np.isfinite(master))

            group.master_filepath = os.path.join(
                self.run_dir, f"Master_{group.group_key}_{group.successful_frames}fr.fits"
            )
            FITSHandler.save_master(
                master, master_err, master_dq,
                group.anchor_filepath, group.master_filepath,
                group.successful_frames, group.meta,
            )

            self._log_telemetry(group, r_data, master, qa_stars, qa_scales, qa_noises)
            del stack, var_stack, cube, varcube, master, master_var, master_err
            gc.collect()

            self.logger.info("  [*] Initiating Astrometry WCS Solve...")
            if not wcs_solver.solve(group, keep_temps=self.keep_temps):
                failed_wcs.append(group)

        if failed_wcs and wcs_solver.global_anchor_ra is not None:
            self.logger.info("\n--- INITIATING WCS RESCUE PASS ---")
            for g in failed_wcs:
                wcs_solver.solve(g, is_rescue=True, keep_temps=self.keep_temps)

        VisualQAGenerator.generate_pdf(self.target_groups, self.run_dir, self.logger, self.config)
        self.logger.info("\n" + "=" * 50)
        self.logger.info("PHASE 2 PIPELINE COMPLETE.")

    def _log_telemetry(self, group, r_data, master, qa_stars, qa_scales, qa_noises):
        """Emit the automated QA telemetry (SNR boost / integration efficiency)."""
        self.logger.info("\n   >>> AUTOMATED QA TELEMETRY <<<")
        if qa_stars:
            self.logger.info(f"   [*] STAR MATCH : Avg: {int(np.mean(qa_stars))} | Min: {int(np.min(qa_stars))}")
        if qa_scales:
            self.logger.info(f"   [*] FLUX SCALE : Avg: {np.mean(qa_scales):.3f}x")

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


def run(input_dir, keep_temps=False, config=None, logger=None, assume_yes=False):
    """Set up and run phase 2 over a directory of calibrated frames."""
    input_dir = os.path.abspath(input_dir)
    if not os.path.exists(input_dir):
        raise SystemExit(f"Directory {input_dir} does not exist.")

    pipeline = IntegrationPipeline(input_dir, keep_temps=keep_temps, config=config, logger=logger)
    pipeline.setup()

    if not assume_yes:
        proceed = input("\nProceed with Unified Integration? [Y/n]: ").strip().lower()
        if proceed not in ("y", ""):
            return
    pipeline.execute()
