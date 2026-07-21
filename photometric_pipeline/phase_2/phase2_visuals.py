# phase2_visuals.py
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.visualization import ZScaleInterval
import warnings

from phase2_math import MathEngine

warnings.filterwarnings('ignore')

class VisualQAGenerator:
    @staticmethod
    def generate_pdf(target_groups, run_dir, logger):
        pdf_path = os.path.join(run_dir, f"Visual_QA_{os.path.basename(run_dir)}.pdf")
        logger.info(f"\n[*] Generating 3-Panel Visual QA PDF: {os.path.basename(pdf_path)}")
        zscale = ZScaleInterval()
        
        valid_groups = [g for g in target_groups if g.master_filepath and g.anchor_filepath]
        if not valid_groups:
            logger.error("[-] No valid Master/Anchor pairs found for PDF generation.")
            return

        with PdfPages(pdf_path) as pdf:
            for group in valid_groups:
                logger.info(f"  -> Rendering: {group.group_key}")
                
                # 1. Load Master Stack
                m_data = fits.getdata(group.master_filepath).astype(np.float32)
                
                # 2. Load Calibrated Anchor
                a_calibrated_raw = fits.getdata(group.anchor_filepath).astype(np.float32)
                calibrated_filename = os.path.basename(group.anchor_filepath)
                
                # 3. Dynamically Find True Raw File
                raw_filename = calibrated_filename.replace("calibrated_", "")
                calibrated_dir = os.path.dirname(group.anchor_filepath)
                parent_dir = os.path.dirname(calibrated_dir)
                raw_dir = os.path.join(parent_dir, "Raw")
                raw_filepath = os.path.join(raw_dir, raw_filename)
                
                if os.path.exists(raw_filepath):
                    a_raw = fits.getdata(raw_filepath).astype(np.float32)
                    raw_title = f"1. True Raw\n({raw_filename})"
                else:
                    a_raw = np.zeros_like(a_calibrated_raw)
                    raw_title = f"1. True Raw\n(FILE NOT FOUND)"
                    logger.warning(f"     -> [Warning] Could not find {raw_filename} in {raw_dir}.")
                
                # Background subtract the single frames for visual brightness matching
                a_raw_bg, _ = MathEngine.extract_2d_background(a_raw)
                a_calibrated_bg, _ = MathEngine.extract_2d_background(a_calibrated_raw)
                
                a_raw_bg += 100.0
                a_calibrated_bg += 100.0
                
                # Get limits from the highest SNR image (the Master)
                vmin, vmax = zscale.get_limits(m_data)
                
                # Create the 1x3 plot
                fig, axes = plt.subplots(1, 3, figsize=(24, 8), dpi=150)
                fig.suptitle(f"Pipeline Progression QA: {group.group_key}", fontsize=18, fontweight='bold', y=0.98)
                
                # Panel 1: Raw
                axes[0].imshow(a_raw_bg, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
                axes[0].set_title(raw_title, fontsize=12)
                axes[0].axis('off')
                
                # Panel 2: Calibrated
                axes[1].imshow(a_calibrated_bg, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
                axes[1].set_title(f"2. Calibrated Frame\n({calibrated_filename})", fontsize=12)
                axes[1].axis('off')
                
                # Panel 3: Master Stack
                im2 = axes[2].imshow(m_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
                axes[2].set_title(f"3. Master Stack\n({group.successful_frames} Frames, WCS: {group.wcs_solved})", fontsize=12)
                axes[2].axis('off')
                
                # Shared Colorbar
                cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
                fig.colorbar(im2, cax=cbar_ax, label='Pixel ADU (Z-Scale Matched)')
                
                plt.subplots_adjust(left=0.03, right=0.9, top=0.88, bottom=0.05, wspace=0.05)
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)