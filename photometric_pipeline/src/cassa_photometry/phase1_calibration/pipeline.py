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
from cassa_photometry.instruments import get_profile
from cassa_photometry.phase1_calibration.data_models import load_standardized_ccds
from cassa_photometry.phase1_calibration.processor import (
    UniversalProcessor, subtract_scaled_dark,
)


def _combine_with_uncertainty(ccd_list, exposure=None, gain=None, bias_subtracted=None):
    """Median-combine frames and attach the standard error of the median.

    The per-pixel uncertainty is ``mad_std`` across the stack divided by
    ``sqrt(N)`` -- an honest, robust estimate of the master-frame noise that
    ``ccdproc`` then propagates into every science frame.

    Provenance the later steps must not have to guess -- the exposure the master
    represents, the detector gain, and whether the bias pedestal has been
    removed -- is recorded in the master's header.
    """
    master = ccdproc.Combiner(ccd_list).median_combine()
    cube = np.array([c.data for c in ccd_list], dtype=np.float64)
    n = cube.shape[0]
    sigma = mad_std(cube, axis=0, ignore_nan=True) / np.sqrt(max(n, 1))
    master.uncertainty = StdDevUncertainty(np.nan_to_num(sigma, nan=np.nanmedian(sigma)))
    master.unit = u.adu
    master.meta["NCOMBINE"] = n
    if exposure is not None:
        master.meta["EXPTIME"] = float(exposure)
    if gain is not None:
        master.meta["GAINVAL"] = float(gain)
    if bias_subtracted is not None:
        master.meta["BIASSUB"] = bool(bias_subtracted)
    return master


def _stack_scatter(master):
    """Per-pixel scatter across the combined stack, recovered from the stored SEM.

    ``_combine_with_uncertainty`` stores ``mad_std / sqrt(N)``; multiplying back
    by ``sqrt(N)`` recovers the frame-to-frame scatter without keeping a second
    full-frame array alive. Returns ``None`` when too few frames went in for the
    scatter to mean anything.
    """
    n = int(master.meta.get("NCOMBINE", 0) or 0)
    if n < 3 or master.uncertainty is None:
        return None
    return master.uncertainty.array * np.sqrt(n)


def _outlier_threshold(values, factor, n_sigma):
    """Upper cut that is both a multiple of the median and a robust sigma clip.

    Taking the larger of the two keeps the test meaningful at both extremes. On a
    detector with appreciable dark current the multiplicative term dominates; on
    a very clean one, where the median sits near zero and a bare multiple would
    flag ordinary scatter across the whole frame, the sigma term takes over.
    """
    median = float(np.nanmedian(values))
    sigma = float(mad_std(values, ignore_nan=True))
    if not np.isfinite(sigma):
        sigma = 0.0
    return max(median * factor, median + n_sigma * sigma)


def build_master_bias(file_list, instrument, config, logger):
    """Median-combine raw bias frames into a master bias (with uncertainty)."""
    if not file_list:
        return None
    ccds = [load_standardized_ccds(f, instrument, config)[0].ccd
            for f in tqdm(file_list, desc="Master Bias", unit="frame")]
    logger.info(f"Median-combining {len(file_list)} bias frames...")
    header = fits.getheader(file_list[0])
    return _combine_with_uncertainty(
        ccds, exposure=instrument.get_exposure(header), gain=instrument.get_gain(header),
    )


def build_master_dark(file_list, instrument, master_bias, config, logger):
    """Build a master dark, bias-subtracting each frame first.

    The master records the exposure time it represents so that
    :func:`~cassa_photometry.phase1_calibration.processor.subtract_scaled_dark`
    can rescale it per frame. A directory holding a *mix* of dark exposures still
    yields one valid master: the frames are normalised to a common exposure
    before combining, which is legitimate once the bias pedestal is gone.
    """
    if not file_list:
        return None
    ccds, exposures = [], []
    for f in tqdm(file_list, desc="Master Dark", unit="frame"):
        std_ccd = load_standardized_ccds(f, instrument, config)[0]
        ccd = std_ccd.ccd
        if master_bias is not None:
            ccd = ccdproc.subtract_bias(ccd, master_bias)
        ccds.append(ccd)
        exposures.append(float(std_ccd.meta["exposure"]))

    ref_exposure = float(np.median(exposures))
    distinct = sorted({round(e, 3) for e in exposures})
    if len(distinct) > 1:
        if master_bias is None:
            logger.warning(
                f"Dark frames have mixed exposures {distinct} but no master bias is "
                f"available, so they cannot be rescaled safely. Using only the "
                f"{ref_exposure:g}s frames."
            )
            kept = [c for c, e in zip(ccds, exposures)
                    if np.isclose(e, ref_exposure, rtol=1e-3)]
            ccds = kept or ccds
        else:
            logger.info(
                f"Dark frames have mixed exposures {distinct}; "
                f"rescaling each to {ref_exposure:g}s before combining."
            )
            for ccd, exposure in zip(ccds, exposures):
                if exposure > 0 and not np.isclose(exposure, ref_exposure, rtol=1e-3):
                    ratio = ref_exposure / exposure
                    ccd.data = ccd.data * ratio
                    if ccd.uncertainty is not None:
                        ccd.uncertainty = StdDevUncertainty(ccd.uncertainty.array * ratio)

    logger.info(f"Median-combining {len(ccds)} dark frames ({ref_exposure:g}s equivalent)...")
    return _combine_with_uncertainty(
        ccds,
        exposure=ref_exposure,
        gain=instrument.get_gain(fits.getheader(file_list[0])),
        bias_subtracted=master_bias is not None,
    )


def build_master_flat(file_list, instrument, master_bias, master_dark, config, logger):
    """Build a master flat: bias- and dark-subtract, normalise *each* frame, combine.

    Dark-subtracting the flats matters little for a few-second exposure on a cold
    detector, but it is not free to skip: a warm detector or a long sky flat
    accumulates real dark current, which would otherwise be normalised into the
    flat and divided back out of every science frame.

    Each frame is then divided by its own median *before* combining, which is the
    order that matters. Twilight sky flats fade as they are taken -- on the
    workshop set the level drifts by 27-57% within a single filter -- and a median
    taken across unnormalised frames is dominated by that drift rather than by the
    per-pixel response. Worse, the stack scatter that becomes the master's
    uncertainty would then be measuring the fading sky, overestimating the flat's
    error by more than an order of magnitude and injecting that inflation into
    every science frame's ERR plane.
    """
    if not file_list:
        return None
    ccds = []
    for f in tqdm(file_list, desc="Master Flat", unit="frame"):
        std_ccd = load_standardized_ccds(f, instrument, config)[0]
        ccd = std_ccd.ccd
        if master_bias is not None:
            ccd = ccdproc.subtract_bias(ccd, master_bias)
        if master_dark is not None:
            ccd = subtract_scaled_dark(ccd, master_dark, std_ccd.meta["exposure"], logger)

        # Convert this frame from ADU into a response map centred on 1.0, so the
        # combine sees only pixel-to-pixel response and not the lamp/sky level.
        level = np.nanmedian(ccd.data)
        if level > 0:
            ccd.data = ccd.data / level
            if ccd.uncertainty is not None:
                ccd.uncertainty = StdDevUncertainty(ccd.uncertainty.array / level)
        else:
            logger.warning(
                f"Flat {os.path.basename(f)} has a non-positive median ({level:.3g}); "
                f"combining it unnormalised."
            )
        ccds.append(ccd)

    levels = [float(np.nanmedian(fits.getdata(f))) for f in file_list]
    logger.info(
        f"Median-combining {len(file_list)} normalised flat frames "
        f"(raw levels {min(levels):.0f}-{max(levels):.0f} ADU)..."
    )
    m_flat = _combine_with_uncertainty(ccds)

    # The combine of already-normalised frames sits at ~1.0; this pins it exactly.
    flat_median = np.nanmedian(m_flat.data)
    if flat_median > 0:
        m_flat.data = m_flat.data / flat_median
        m_flat.uncertainty = StdDevUncertainty(m_flat.uncertainty.array / flat_median)
    return m_flat


def _build_bpm(master_flats, master_dark, master_bias, config, logger):
    """Derive a bad-pixel mask (True == bad) from every master available.

    Each calibration product sees a defect class the others are blind to, so all
    three are consulted and the results OR-ed together:

    * **flats** -- anomalous *sensitivity*: dead, low-QE and dust-shadowed pixels.
    * **darks** -- *hot* pixels, i.e. anomalous dark current. A short flat cannot
      see these at all: a few seconds of even severe dark current is a rounding
      error against the lamp signal.
    * **biases** -- *unstable* pixels whose frame-to-frame scatter far exceeds the
      rest of the detector (random-telegraph / flickering pixels).

    Thresholds are expressed relative to each frame's own median rather than in
    absolute ADU or e-/s, so they carry across detectors without retuning. Any
    test whose master is missing is skipped, and a factor of 0 disables one.
    """
    cfg = config.phase1
    bpm = None
    reasons = {}

    def _add(mask, label):
        nonlocal bpm
        if mask is None:
            return
        mask = np.asarray(mask, dtype=bool)
        reasons[label] = int(mask.sum())
        bpm = mask if bpm is None else (bpm | mask)

    # 1. Sensitivity defects, from the normalised master flats.
    flat_bad = None
    for m_flat in master_flats.values():
        bad = (
            (m_flat.data < cfg.bpm_flat_low)
            | (m_flat.data > cfg.bpm_flat_high)
            | ~np.isfinite(m_flat.data)
        )
        flat_bad = bad if flat_bad is None else (flat_bad | bad)
    _add(flat_bad, "dead/low-QE (flat)")

    # 2. Hot pixels, from the master dark's implied dark-current rate.
    if master_dark is not None and cfg.bpm_dark_rate_factor:
        exposure = float(master_dark.meta.get("EXPTIME") or 0.0)
        gain = float(master_dark.meta.get("GAINVAL") or 1.0)
        if exposure > 0:
            rate = master_dark.data * gain / exposure           # e-/s
            threshold = _outlier_threshold(rate, cfg.bpm_dark_rate_factor, cfg.bpm_dark_sigma)
            logger.info(
                f"Dark current: median {np.nanmedian(rate):.4f} e-/s, "
                f"hot-pixel threshold {threshold:.4f} e-/s."
            )
            _add(rate > threshold, "hot (dark)")
        else:
            logger.warning("Master dark has no exposure time; skipping the hot-pixel test.")

    # 3. Unstable pixels, from the scatter across the bias stack.
    if master_bias is not None and cfg.bpm_bias_noise_factor:
        scatter = _stack_scatter(master_bias)
        if scatter is None:
            logger.info("Too few bias frames to measure pixel stability; skipping that test.")
        elif np.nanmedian(scatter) > 0:
            threshold = _outlier_threshold(
                scatter, cfg.bpm_bias_noise_factor, cfg.bpm_bias_noise_sigma
            )
            logger.info(
                f"Bias stability: median scatter {np.nanmedian(scatter):.3f} ADU, "
                f"unstable-pixel threshold {threshold:.3f} ADU."
            )
            _add(scatter > threshold, "unstable (bias)")

    if bpm is not None:
        detail = ", ".join(f"{n} {label}" for label, n in reasons.items())
        logger.info(f"Bad-pixel mask: {int(np.sum(bpm))} pixels flagged ({detail}).")
    return bpm


def run(data_dir, output_dir, config=None, logger=None, instrument=None):
    """Run phase 1 (ISR) over a directory of raw FITS frames.

    ``instrument`` defaults to the profile named by ``config.instrument``.
    """
    config = config or load_config()
    logger = logger or get_logger("cassa_calibrate")
    instrument = instrument or get_profile(config.instrument)

    logger.info(f"Initializing pipeline: {instrument.name}")
    logger.info(f"Input:  {data_dir}")
    logger.info(f"Output: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    all_files = glob.glob(os.path.join(data_dir, "*.f*t*"))
    if not all_files:
        logger.error(f"No FITS files found in {data_dir}")
        return

    categorized = {"bias": [], "dark": [], "flat": {}, "science": {}}
    precalibrated = []
    for f in tqdm(all_files, desc="Scanning headers", unit="file"):
        try:
            with fits.open(f) as hdul:
                header = hdul[0].header
                # Reducing an already-reduced frame produces a plausible-looking
                # image with a wrong error budget, and nothing downstream can
                # tell. Catch it here, once, before any work is done.
                reason = instrument.already_calibrated(header)
                if reason:
                    precalibrated.append((os.path.basename(f), reason))
                    if not config.phase1.allow_precalibrated:
                        continue
                img_type = instrument.get_image_type(header)
                if img_type in ("flat", "science"):
                    filt = instrument.get_filter(header)
                    categorized[img_type].setdefault(filt, []).append(f)
                else:
                    categorized[img_type].append(f)
        except Exception as exc:
            logger.warning(f"Error reading {os.path.basename(f)}: {exc}")

    if precalibrated:
        verb = "Reducing" if config.phase1.allow_precalibrated else "Skipping"
        logger.warning(
            f"{verb} {len(precalibrated)} frame(s) that report prior calibration: "
            + "; ".join(f"{name} ({why})" for name, why in precalibrated[:5])
            + (" ..." if len(precalibrated) > 5 else "")
        )
        if not config.phase1.allow_precalibrated:
            logger.warning(
                "Raw frames should have an empty CALSTAT and be in ADU. Set "
                "phase1.allow_precalibrated: true to reduce them anyway."
            )

    logger.info("Building master frames...")
    m_bias = build_master_bias(categorized["bias"], instrument, config, logger)
    m_dark = build_master_dark(categorized["dark"], instrument, m_bias, config, logger)

    m_flats = {}
    for filt, files in categorized["flat"].items():
        m_flats[filt] = build_master_flat(files, instrument, m_bias, m_dark, config, logger)

    # Flat reuse is a property of the instrument, so the profile declares it.
    for missing, proxy in instrument.flat_proxies().items():
        if missing not in m_flats and proxy in m_flats:
            logger.info(f"Missing {missing} flat; using the {proxy} master flat as proxy.")
            m_flats[missing] = m_flats[proxy]

    bpm = _build_bpm(m_flats, m_dark, m_bias, config, logger)

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
