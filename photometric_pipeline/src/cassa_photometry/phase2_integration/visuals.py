"""Three-panel visual QA PDF (raw -> calibrated -> master stack) for phase 2."""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from matplotlib.backends.backend_pdf import PdfPages

from cassa_photometry.config import load_config
from cassa_photometry.paths import find_raw_frames
from cassa_photometry.phase2_integration.math_utils import MathEngine

# Scoped, not global. A module-level filterwarnings("ignore") silences warnings
# for the whole *process* -- including astropy's WCS and FITS-verification
# warnings, which are exactly the ones a reduction wants to see.
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")


def _raw_candidates(anchor_path, calibrated_filename, header, configured_dir):
    """Yield ``(filename, directory)`` pairs to try, best provenance first.

    The raw tree is a CLI argument and routinely sits outside the work
    directory, so guessing a sibling of the calibrated directory only works for
    one particular layout. Phase 1 stamps ``RAWFILE``/``RAWDIR`` into every
    calibrated header precisely so this does not have to be guessed; the
    directory guesses below are the fallback for frames calibrated before those
    cards existed, or for a raw tree that has since moved.
    """
    stamped = header.get("RAWFILE")
    name = str(stamped).strip() if stamped else calibrated_filename.replace("calibrated_", "", 1)

    calibrated_dir = os.path.dirname(os.path.abspath(anchor_path))
    work_dir = os.path.dirname(calibrated_dir)

    dirs = []
    if configured_dir:
        dirs.append(os.path.abspath(os.path.expanduser(configured_dir)))
    stamped_dir = header.get("RAWDIR")
    if stamped_dir:
        dirs.append(str(stamped_dir).strip())
    # Conventional layouts, in decreasing confidence. Both cases of "raw" are
    # tried because the directory name comes from the user, not the pipeline.
    for base in (work_dir, os.path.dirname(work_dir)):
        for leaf in ("Raw", "raw", "RAW"):
            dirs.append(os.path.join(base, leaf))

    seen = set()
    for directory in dirs:
        if directory and directory not in seen:
            seen.add(directory)
            yield name, directory


def find_raw_frame(anchor_path, calibrated_filename, header, configured_dir=None):
    """Locate the raw frame behind a calibrated one, or ``(None, searched)``.

    Each candidate directory is tried as the frame's own directory first and
    then as the root of a tree, because the acquisition software files a night
    under ``<date>/<TYPE>/<target>`` and a raw argument naming the night's root
    is the normal case rather than the exception.
    """
    searched = []
    for name, directory in _raw_candidates(anchor_path, calibrated_filename,
                                           header, configured_dir):
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            return candidate, searched
        nested = _find_below(directory, name)
        if nested:
            return nested, searched
        searched.append(directory)
    return None, searched


def _find_below(directory, name):
    """First frame called ``name`` anywhere under ``directory``, or ``None``."""
    if not os.path.isdir(directory):
        return None
    for path in find_raw_frames(directory):
        if os.path.basename(path) == name:
            return path
    return None


class VisualQAGenerator:
    @staticmethod
    def generate_pdf(target_groups, run_dir, logger, config=None):
        config = config or load_config()
        box = config.phase2.background_box
        filt = config.phase2.background_filter

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

                # fits.getdata returns the primary (SCI) plane of the MEF files.
                m_data = fits.getdata(group.master_filepath).astype(np.float32)
                with fits.open(group.anchor_filepath) as hdul:
                    a_calibrated_raw = np.asarray(hdul[0].data, dtype=np.float32)
                    anchor_header = hdul[0].header.copy()
                calibrated_filename = os.path.basename(group.anchor_filepath)

                raw_filepath, searched = find_raw_frame(
                    group.anchor_filepath, calibrated_filename, anchor_header,
                    config.phase2.raw_dir,
                )

                if raw_filepath:
                    a_raw = fits.getdata(raw_filepath).astype(np.float32)
                    raw_title = f"1. True Raw\n({os.path.basename(raw_filepath)})"
                else:
                    a_raw = np.zeros_like(a_calibrated_raw)
                    raw_title = "1. True Raw\n(FILE NOT FOUND)"
                    raw_name = str(anchor_header.get(
                        "RAWFILE", calibrated_filename.replace("calibrated_", "", 1))).strip()
                    logger.warning(
                        f"     -> [Warning] Could not find {raw_name}; searched "
                        f"{', '.join(searched)}. Set phase2.raw_dir to the raw "
                        f"directory (frames calibrated before RAWDIR was stamped "
                        f"carry no provenance)."
                    )

                a_raw_bg, _ = MathEngine.extract_2d_background(a_raw, box, filt)
                a_calibrated_bg, _ = MathEngine.extract_2d_background(a_calibrated_raw, box, filt)

                vmin, vmax = zscale.get_limits(m_data)

                fig, axes = plt.subplots(1, 3, figsize=(24, 8), dpi=150)
                fig.suptitle(f"Pipeline Progression QA: {group.group_key}",
                            fontsize=18, fontweight="bold", y=0.98)

                axes[0].imshow(a_raw_bg, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
                axes[0].set_title(raw_title, fontsize=12)
                axes[0].axis("off")

                axes[1].imshow(a_calibrated_bg, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
                axes[1].set_title(f"2. Calibrated Frame\n({calibrated_filename})", fontsize=12)
                axes[1].axis("off")

                im2 = axes[2].imshow(m_data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
                axes[2].set_title(
                    f"3. Master Stack\n({group.successful_frames} Frames, WCS: {group.wcs_solved})",
                    fontsize=12,
                )
                axes[2].axis("off")

                cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
                fig.colorbar(im2, cax=cbar_ax, label="Pixel electrons (Z-Scale Matched)")

                plt.subplots_adjust(left=0.03, right=0.9, top=0.88, bottom=0.05, wspace=0.05)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
