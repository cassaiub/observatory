"""Phase 1 orchestrator: build master calibrations and reduce science frames.

Produces multi-extension ``calibrated_*.fits`` files (SCI/ERR/DQ, electrons) in
the output directory, conventionally ``<work>/phase1``. Master frames are built
*with* uncertainty (the standard error of the robust combine) so the science
error budget includes the calibration noise. The masters themselves are used and
discarded -- only the calibrated science frames are written.
"""

import os
import glob

import numpy as np
import ccdproc
import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.stats import mad_std
from tqdm import tqdm

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.fits_utils import write_mef
from cassa_photometry.instruments import ITelescopeNetworkProfile
from cassa_photometry.phase1_calibration.data_models import load_standardized_ccds
from cassa_photometry.phase1_calibration.processor import UniversalProcessor


def _combine_with_uncertainty(ccd_list):
    """Median-combine frames and attach the standard error of the median.

    The per-pixel uncertainty is ``mad_std`` across the stack divided by
    ``sqrt(N)`` -- an honest, robust estimate of the master-frame noise that
    ``ccdproc`` then propagates into every science frame.
    """
    master = ccdproc.Combiner(ccd_list).median_combine()
    cube = np.array([c.data for c in ccd_list], dtype=np.float64)
    n = cube.shape[0]
    sigma = mad_std(cube, axis=0, ignore_nan=True) / np.sqrt(max(n, 1))
    master.uncertainty = StdDevUncertainty(np.nan_to_num(sigma, nan=np.nanmedian(sigma)))
    master.unit = u.adu
    return master


def build_master_bias(file_list, instrument, config, logger):
    """Median-combine raw bias frames into a master bias (with uncertainty)."""
    if not file_list:
        return None
    ccds = [load_standardized_ccds(f, instrument, config)[0].ccd
            for f in tqdm(file_list, desc="Master Bias", unit="frame")]
    logger.info(f"Median-combining {len(file_list)} bias frames...")
    return _combine_with_uncertainty(ccds)


def build_master_dark(file_list, instrument, master_bias, config, logger):
    """Build a master dark, bias-subtracting each frame first."""
    if not file_list:
        return None
    ccds = []
    for f in tqdm(file_list, desc="Master Dark", unit="frame"):
        ccd = load_standardized_ccds(f, instrument, config)[0].ccd
        if master_bias is not None:
            ccd = ccdproc.subtract_bias(ccd, master_bias)
        ccds.append(ccd)
    logger.info(f"Median-combining {len(file_list)} dark frames...")
    return _combine_with_uncertainty(ccds)


def build_master_flat(file_list, instrument, master_bias, config, logger):
    """Build a master flat: bias-subtract, combine, and normalise to 1.0."""
    if not file_list:
        return None
    ccds = []
    for f in tqdm(file_list, desc="Master Flat", unit="frame"):
        ccd = load_standardized_ccds(f, instrument, config)[0].ccd
        if master_bias is not None:
            ccd = ccdproc.subtract_bias(ccd, master_bias)
        ccds.append(ccd)
    logger.info(f"Median-combining {len(file_list)} flat frames...")
    m_flat = _combine_with_uncertainty(ccds)

    flat_median = np.nanmedian(m_flat.data)
    if flat_median > 0:
        logger.info(f"Normalising master flat (median ADU: {flat_median:.1f})")
        m_flat.data = m_flat.data / flat_median
        m_flat.uncertainty = StdDevUncertainty(m_flat.uncertainty.array / flat_median)
    return m_flat


def _build_bpm(master_flats, config, logger):
    """Derive a bad-pixel mask (True == bad) from the normalised master flats."""
    cfg = config.phase1
    bpm = None
    for m_flat in master_flats.values():
        bad = (
            (m_flat.data < cfg.bpm_flat_low)
            | (m_flat.data > cfg.bpm_flat_high)
            | ~np.isfinite(m_flat.data)
        )
        bpm = bad if bpm is None else (bpm | bad)
    if bpm is not None:
        logger.info(f"Bad-pixel mask: {int(np.sum(bpm))} pixels flagged.")
    return bpm


def run(data_dir, output_dir, config=None, logger=None):
    """Run phase 1 (ISR) over a directory of raw FITS frames."""
    config = config or load_config()
    logger = logger or get_logger("cassa_calibrate")
    instrument = ITelescopeNetworkProfile()

    logger.info(f"Initializing pipeline: {instrument.name}")
    logger.info(f"Input:  {data_dir}")
    logger.info(f"Output: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    all_files = glob.glob(os.path.join(data_dir, "*.f*t*"))
    if not all_files:
        logger.error(f"No FITS files found in {data_dir}")
        return

    categorized = {"bias": [], "dark": [], "flat": {}, "science": {}}
    for f in tqdm(all_files, desc="Scanning headers", unit="file"):
        try:
            with fits.open(f) as hdul:
                header = hdul[0].header
                img_type = instrument.get_image_type(header)
                if img_type in ("flat", "science"):
                    filt = instrument.get_filter(header)
                    categorized[img_type].setdefault(filt, []).append(f)
                else:
                    categorized[img_type].append(f)
        except Exception as exc:
            logger.warning(f"Error reading {os.path.basename(f)}: {exc}")

    logger.info("Building master frames...")
    m_bias = build_master_bias(categorized["bias"], instrument, config, logger)
    m_dark = build_master_dark(categorized["dark"], instrument, m_bias, config, logger)

    m_flats = {}
    for filt, files in categorized["flat"].items():
        m_flats[filt] = build_master_flat(files, instrument, m_bias, config, logger)

    if "RED" not in m_flats and "LUMINANCE" in m_flats:
        logger.info("Missing RED flat; mapping LUMINANCE master flat as proxy for RED.")
        m_flats["RED"] = m_flats["LUMINANCE"]

    bpm = _build_bpm(m_flats, config, logger)

    processor = UniversalProcessor(
        master_bias=m_bias, master_dark=m_dark, master_flats=m_flats,
        bpm=bpm, config=config, logger=logger,
    )

    if not categorized["science"]:
        logger.warning("No science frames found to process.")
        return

    for filt, sci_files in categorized["science"].items():
        logger.info(f"Processing {filt} filter frames ({len(sci_files)})")
        for sci_file in tqdm(sci_files, desc=f"Calibrating {filt}", unit="img"):
            standard_ccds = load_standardized_ccds(sci_file, instrument, config)
            for frame in processor.process_science_frame(standard_ccds):
                out_name = f"calibrated_{os.path.basename(sci_file)}"
                out_path = os.path.join(output_dir, out_name)
                write_mef(
                    out_path,
                    sci=frame.ccd.data,
                    err=frame.ccd.uncertainty.array if frame.ccd.uncertainty is not None else None,
                    dq=frame.dq,
                    header=frame.ccd.header,
                    history="PHASE 1: ISR (bias/dark/flat, CR, gain); SCI/ERR/DQ in electrons",
                )

    logger.info("Phase 1 complete. Calibrated SCI/ERR/DQ frames written to output.")
