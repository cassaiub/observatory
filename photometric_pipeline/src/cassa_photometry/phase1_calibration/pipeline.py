"""Phase 1 orchestrator: build master calibrations and reduce science frames.

Produces multi-extension ``calibrated_*.fits`` files (SCI/ERR/DQ, electrons) in
the output directory, conventionally ``<work>/phase1``. Master frames are built
*with* uncertainty (the standard error of the robust combine) so the science
error budget includes the calibration noise. The masters themselves are used and
discarded -- only the calibrated science frames are written.
"""

import os

import astropy.units as u
import ccdproc
import numpy as np
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.stats import mad_std
from tqdm import tqdm

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import CALVERS, open_fits, write_mef
from cassa_photometry.instruments import get_profile
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import find_raw_frames, raw_tree_summary
from cassa_photometry.phase1_calibration.data_models import load_standardized_ccds
from cassa_photometry.phase1_calibration.processor import (
    UniversalProcessor,
    subtract_scaled_dark,
)

#: Standard error of the *median* is wider than that of the mean by this factor
#: for Gaussian noise (sqrt(pi/2)). Ignoring it understates every master's
#: uncertainty by 25%, which then propagates into every science frame's ERR.
MEDIAN_SEM_FACTOR = np.sqrt(np.pi / 2.0)

#: Below this many frames the frame-to-frame scatter is too poorly determined to
#: use as an uncertainty: with N=1 it is identically zero, which would give the
#: master an ERR of 0 and let it contribute nothing at all to the error budget.
MIN_FRAMES_FOR_EMPIRICAL_SIGMA = 3


def _bias_level(master_bias):
    """Median bias pedestal in ADU, or 0.0 when there is no master bias.

    Handed to the frame loader so the Poisson term is computed from collected
    charge rather than from charge plus an electronic offset.
    """
    if master_bias is None:
        return 0.0
    import numpy as _np

    level = float(_np.nanmedian(_np.asarray(master_bias.data, dtype=float)))
    return level if _np.isfinite(level) else 0.0


def _combine_with_uncertainty(ccd_list, exposure=None, gain=None, bias_subtracted=None,
                              sigma_clip=3.0, logger=None):
    """Combine frames with outlier rejection and attach an honest uncertainty.

    Two things here are easy to get subtly wrong, and both were:

    * **The uncertainty is the standard error of the median, not of the mean.**
      ``mad_std / sqrt(N)`` is the SEM of the *mean*; the median's is wider by
      ``sqrt(pi/2)`` for Gaussian noise. Understating it by 25% understates every
      science frame's ERR that this master touches.
    * **A stack of one or two frames has no usable scatter.** ``mad_std`` over an
      axis of length 1 is exactly zero, so a master built from a single bias
      would claim perfect knowledge and contribute nothing to the error budget.
      Below :data:`MIN_FRAMES_FOR_EMPIRICAL_SIGMA` the propagated per-frame
      uncertainty is used instead, and the situation is logged.

    Provenance the later steps must not have to guess -- the exposure the master
    represents, the detector gain, and whether the bias pedestal has been
    removed -- is recorded in the master's header.
    """
    combiner = ccdproc.Combiner(ccd_list)
    n = len(ccd_list)

    # Reject cosmic rays and satellite trails from the calibration stack itself.
    # A median is robust to a few outliers but not immune, and nothing else in
    # phase 1 looks at calibration frames for them.
    n_rejected = 0
    if sigma_clip and n >= MIN_FRAMES_FOR_EMPIRICAL_SIGMA:
        combiner.sigma_clipping(low_thresh=sigma_clip, high_thresh=sigma_clip,
                                func=np.ma.median, dev_func=mad_std)
        n_rejected = int(np.count_nonzero(combiner.data_arr.mask))

    master = combiner.median_combine()
    cube = np.array([c.data for c in ccd_list], dtype=np.float64)

    if n >= MIN_FRAMES_FOR_EMPIRICAL_SIGMA:
        sigma = MEDIAN_SEM_FACTOR * mad_std(cube, axis=0, ignore_nan=True) / np.sqrt(n)
        sigma = np.nan_to_num(sigma, nan=float(np.nanmedian(sigma)))
    else:
        sigma = _propagated_sigma(ccd_list, cube, n, logger)

    master.uncertainty = StdDevUncertainty(sigma)
    master.unit = u.adu
    master.meta["NCOMBINE"] = n
    master.meta["NREJECT"] = n_rejected
    if exposure is not None:
        master.meta["EXPTIME"] = float(exposure)
    if gain is not None:
        master.meta["GAINVAL"] = float(gain)
    if bias_subtracted is not None:
        master.meta["BIASSUB"] = bool(bias_subtracted)
    if logger is not None and n_rejected:
        logger.info(
            "    %d pixel value(s) sigma-clipped from the %d-frame stack (%.3f%%).",
            n_rejected, n, 100.0 * n_rejected / cube.size,
        )
    return master


def _propagated_sigma(ccd_list, cube, n, logger=None):
    """Uncertainty for a stack too small to measure its own scatter."""
    if logger is not None:
        logger.warning(
            "    Only %d frame(s) in this master: the frame-to-frame scatter is not "
            "measurable, so the propagated per-frame uncertainty is used instead. "
            "Take at least %d calibration frames.",
            n, MIN_FRAMES_FOR_EMPIRICAL_SIGMA,
        )
    variances = [
        np.asarray(c.uncertainty.array, dtype=float) ** 2
        for c in ccd_list if c.uncertainty is not None
    ]
    if variances:
        # Median of N frames: SEM of the mean, widened for the median.
        return MEDIAN_SEM_FACTOR * np.sqrt(np.sum(variances, axis=0)) / n
    # Nothing to propagate: fall back to the whole-frame scatter, which at least
    # is not zero.
    return np.full(cube.shape[1:], float(mad_std(cube, ignore_nan=True)) or 1.0)


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
        sigma_clip=config.phase1.master_sigma_clip, logger=logger,
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
        std_ccd = load_standardized_ccds(
            f, instrument, config, bias_level=_bias_level(master_bias)
        )[0]
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
            kept = [c for c, e in zip(ccds, exposures, strict=True)
                    if np.isclose(e, ref_exposure, rtol=1e-3)]
            ccds = kept or ccds
        else:
            logger.info(
                f"Dark frames have mixed exposures {distinct}; "
                f"rescaling each to {ref_exposure:g}s before combining."
            )
            for ccd, exposure in zip(ccds, exposures, strict=True):
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
        sigma_clip=config.phase1.master_sigma_clip,
        logger=logger,
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
        std_ccd = load_standardized_ccds(
            f, instrument, config, bias_level=_bias_level(master_bias)
        )[0]
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
    m_flat = _combine_with_uncertainty(
        ccds, sigma_clip=config.phase1.master_sigma_clip, logger=logger
    )

    # The combine of already-normalised frames sits at ~1.0; this pins it exactly.
    flat_median = np.nanmedian(m_flat.data)
    if flat_median > 0:
        m_flat.data = m_flat.data / flat_median
        m_flat.uncertainty = StdDevUncertainty(m_flat.uncertainty.array / flat_median)
    return m_flat


def _flat_relative_response(data, smooth_px):
    """A flat divided by a smoothed copy of itself: pixel response alone.

    Removes vignetting and illumination gradients, which are optics rather than
    detector defects, so a bad-pixel threshold measures what it claims to.
    Falls back to the raw flat when smoothing is disabled or unavailable.
    """
    data = np.asarray(data, dtype=float)
    if not smooth_px or smooth_px <= 0:
        return data
    try:
        from scipy.ndimage import median_filter

        size = max(int(smooth_px), 3)
        smooth = median_filter(np.nan_to_num(data, nan=1.0), size=size, mode="nearest")
    except Exception:
        return data
    with np.errstate(invalid="ignore", divide="ignore"):
        relative = np.where(smooth > 0, data / smooth, np.nan)
    return relative


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
    #
    # Measured against the flat's own *large-scale* level, not against 1.0.
    # Vignetting and an off-axis illumination gradient are properties of the
    # optics, not defects: on a strongly vignetted train the corners can sit at
    # 0.5-0.6 of the centre, and a flat threshold would condemn them wholesale.
    # Dividing by a smoothed copy leaves only the pixel-to-pixel component,
    # which is what "is this pixel bad" actually means.
    flat_bad = None
    for m_flat in master_flats.values():
        relative = _flat_relative_response(m_flat.data, cfg.bpm_flat_smooth_px)
        bad = (
            (relative < cfg.bpm_flat_low)
            | (relative > cfg.bpm_flat_high)
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


def _record_frame_fwhm(header, data, config, instrument):
    """Measure and stamp this frame's delivered PSF.

    ``FWHMPX`` is what the pipeline measured from the stars in *this* image: the
    delivered point-spread function, in this band, at this airmass, including
    guiding and optics. It is not ``SEEING``, which acquisition writes from the
    DIMM and which is atmospheric, zenith-corrected and at a reference
    wavelength -- using one as the other under-sizes photometric apertures by
    roughly 30%.

    Two things about the number, both measured against simulated truth:

    * It is a **field average**, not the on-axis PSF. The stars it fits are
      spread across the frame, and an 8-inch Newtonian's PSF broadens off-axis;
      on the simulated CASSA field that makes ``FWHMPX`` 14% larger than the
      centre value. That is the honest quantity for sizing apertures over a whole
      frame, but it is why a *scalar* aperture correction is wrong in the corners.
    * A Gaussian fitted to a Moffat core reads about **5% high**, because the
      wings inside the fitting box pull the Gaussian wider. Apertures sized from
      it are therefore slightly generous, which is the safe direction.

    A star-poor frame simply gets no card, and every consumer falls back.
    """
    from cassa_photometry.phase4_diagnostics.psf import estimate_fwhm

    try:
        pixel_scale = instrument.get_pixel_scale(header)
        result = estimate_fwhm(
            np.asarray(data, dtype=float),
            fwhm_guess=config.phase4.fwhm_guess,
            threshold=config.phase4.detection_threshold,
            pixscale=pixel_scale,
            max_stars=config.phase4.max_stars,
            cutout=config.phase4.cutout,
        )
    except Exception:
        return

    fwhm = result.get("fwhm_px") if result else None
    if fwhm is None or not np.isfinite(fwhm) or fwhm <= 0:
        return
    header["FWHMPX"] = (float(fwhm), "[pixel] Measured delivered PSF FWHM")
    if result.get("fwhm_arcsec"):
        header["FWHMASEC"] = (float(result["fwhm_arcsec"]), "[arcsec] Measured PSF FWHM")
    if result.get("n_used"):
        header["FWHMNSTR"] = (int(result["n_used"]), "Stars used for FWHMPX")
    if result.get("ellipticity") is not None and np.isfinite(result["ellipticity"]):
        header["FWHMELL"] = (float(result["ellipticity"]), "Median stellar ellipticity")


def _setup_of(header, path):
    """The detector configuration a frame was taken under."""
    return {
        "path": path,
        "binning": (header.get("XBINNING"), header.get("YBINNING")),
        "shape": (header.get("NAXIS1"), header.get("NAXIS2")),
        "temperature": header.get("CCD-TEMP", header.get("SET-TEMP")),
        "date": str(header.get("DATE-OBS", ""))[:10],
        "flat_type": str(header.get("FLATTYPE", "") or "").strip().lower(),
    }


def _report_setup_mismatches(setups, config, logger):
    """Warn where calibrations were taken under different conditions than the science.

    Masters are built from every frame of a type in the directory, whatever
    binning, sensor temperature or night it came from. That is usually what the
    user intends, and occasionally it is a silent error -- a dark stack at a
    different set point removes the wrong thermal signal, and a flat at a
    different binning cannot be applied at all. Reporting it costs nothing and
    is the difference between a puzzling result and an obvious one.
    """
    science = setups.get("science") or []
    if not science:
        return

    def _distinct(entries, key):
        return {e[key] for e in entries if e[key] not in (None, "", (None, None))}

    science_binning = _distinct(science, "binning")
    science_temps = _distinct(science, "temperature")

    for kind in ("bias", "dark", "flat"):
        entries = setups.get(kind) or []
        if not entries:
            continue

        binning = _distinct(entries, "binning")
        if science_binning and binning and not (binning & science_binning):
            logger.warning(
                "Master %s frames are binned %s but the science frames are %s. "
                "They cannot be applied and the frames will be skipped.",
                kind, sorted(binning), sorted(science_binning),
            )

        temps = _distinct(entries, "temperature")
        if science_temps and temps:
            try:
                spread = max(abs(float(t) - float(s))
                             for t in temps for s in science_temps)
                if spread > config.phase1.calibration_temp_tolerance_c:
                    logger.warning(
                        "Master %s frames were taken up to %.1f C from the science "
                        "frames' sensor temperature. Dark current roughly doubles "
                        "every 6-7 C, so a mismatched dark removes the wrong "
                        "thermal signal.", kind, spread,
                    )
            except (TypeError, ValueError):
                pass

        if len(_distinct(entries, "date")) > 1:
            logger.info(
                "Master %s frames span %d different dates; they will be combined "
                "into one master.", kind, len(_distinct(entries, "date")),
            )

    # Sky flats correct illumination as well as pixel response; dome and panel
    # flats correct response only, leaving a large-scale gradient that becomes a
    # position-dependent zero point downstream.
    flat_types = {e["flat_type"] for e in (setups.get("flat") or []) if e["flat_type"]}
    if flat_types and not flat_types <= {"sky", "twilight"}:
        logger.warning(
            "FLATTYPE reports %s flats. These correct pixel response but not "
            "illumination, so a large-scale gradient survives into the science "
            "frames and becomes a position-dependent zero point. No illumination "
            "correction is implemented; sky flats are the safer choice.",
            "/".join(sorted(flat_types)),
        )
    elif not flat_types and setups.get("flat"):
        logger.info(
            "Master flats carry no FLATTYPE card, so the pipeline cannot tell a "
            "sky flat from a dome or panel flat. Write FLATTYPE at acquisition."
        )


def _calibrated_names(science_files, data_dir, logger):
    """``{raw path: output filename}``, unique within one run.

    Output is one flat directory, so the name is the raw basename with a
    ``calibrated_`` prefix. A nested raw tree can repeat a basename across
    directories -- two nights of the same target, or the same frame filed under
    two targets -- and the second write would then silently replace the first.
    Collisions are disambiguated by the raw subdirectory instead, and reported,
    because losing frames without a word is the worse failure.
    """
    names, taken, collided = {}, {}, False
    for path in science_files:
        name = f"calibrated_{os.path.basename(path)}"
        if name in taken:
            collided = True
            for candidate in (path, taken[name]):
                relative = os.path.relpath(os.path.dirname(os.path.abspath(candidate)),
                                           os.path.abspath(data_dir))
                prefix = relative.replace(os.sep, "_").strip("._") or "root"
                names[candidate] = f"calibrated_{prefix}_{os.path.basename(candidate)}"
        else:
            taken[name] = path
            names[path] = name
    if collided and logger is not None:
        logger.warning(
            "Raw frames share a filename across directories; those outputs carry "
            "their raw subdirectory in the name so none is overwritten."
        )
    return names


def run(data_dir, output_dir, config=None, logger=None, instrument=None):
    """Run phase 1 (ISR) over a directory of raw FITS frames.

    ``instrument`` defaults to the profile named by ``config.instrument``.
    """
    config = config or load_config()
    own_logger = logger is None
    logger = logger or get_logger("cassa_calibrate")
    instrument = instrument or get_profile(config.instrument, config=config)

    logger.info(f"Initializing pipeline: {instrument.name}")
    logger.info(f"Input:  {data_dir}")
    logger.info(f"Output: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    if own_logger:
        # Phase 2 has always kept a run log; phases 1/3/4 only logged to the
        # console, so a batch run left no record of what it did. get_logger is
        # idempotent, so this just adds the file handler.
        logger = get_logger("cassa_calibrate", run_dir=output_dir)

    # The raw tree is whatever the acquisition software wrote -- a flat
    # directory, or a night sorted into <date>/<TYPE>/<target> -- so it is
    # walked rather than listed. What each frame *is* still comes from its
    # header, not from the folder it was filed under.
    all_files = find_raw_frames(data_dir)
    if not all_files:
        logger.error(f"No FITS files found in {data_dir} or below it")
        return

    layout = raw_tree_summary(all_files, data_dir)
    if list(layout) != ["."]:
        logger.info("Raw tree holds %d frame(s) in %d directories:", len(all_files), len(layout))
        for relative, count in layout.items():
            logger.info("    %-40s %4d", relative, count)

    categorized = {"bias": [], "dark": [], "flat": {}, "science": {}}
    precalibrated = []
    #: Detector setup per frame, so calibrations taken under different
    #: conditions than the science frames are reported rather than applied.
    setups = {}
    for f in tqdm(all_files, desc="Scanning headers", unit="file"):
        try:
            with open_fits(f) as hdul:
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
                setups.setdefault(img_type, []).append(_setup_of(header, f))
                if img_type in ("flat", "science"):
                    filt = instrument.get_filter(header)
                    categorized[img_type].setdefault(filt, []).append(f)
                else:
                    categorized[img_type].append(f)
        except Exception as exc:
            logger.warning(f"Error reading {os.path.basename(f)}: {exc}")

    _report_setup_mismatches(setups, config, logger)

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

    skipped = config.phase1.steps.skipped()
    if skipped:
        logger.warning(
            "Deliberately skipping: %s. The calibrated frames record this as "
            "CALSKIP so a partially-reduced frame cannot pass for a finished one.",
            ", ".join(skipped),
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

    bpm = (_build_bpm(m_flats, m_dark, m_bias, config, logger)
           if config.phase1.steps.enabled("bad_pixel_mask") else None)

    processor = UniversalProcessor(
        master_bias=m_bias, master_dark=m_dark, master_flats=m_flats,
        bpm=bpm, config=config, logger=logger,
    )

    if not categorized["science"]:
        logger.warning("No science frames found to process.")
        return

    # The pedestal the science frames' error budget is measured against.
    science_bias_level = _bias_level(m_bias)

    out_names = _calibrated_names(
        [f for files in categorized["science"].values() for f in files], data_dir, logger)

    for filt, sci_files in categorized["science"].items():
        logger.info(f"Processing {filt} filter frames ({len(sci_files)})")
        for sci_file in tqdm(sci_files, desc=f"Calibrating {filt}", unit="img"):
            standard_ccds = load_standardized_ccds(
                sci_file, instrument, config, bias_level=science_bias_level
            )
            for frame in processor.process_science_frame(standard_ccds):
                out_path = os.path.join(output_dir, out_names[sci_file])
                # Provenance: which raw frame this came from, and where it lived.
                # Later phases need the original to display or re-read it, and
                # the raw directory is not recoverable from the output tree --
                # it is a CLI argument that may sit anywhere on disk.
                frame.ccd.header["RAWFILE"] = (
                    os.path.basename(sci_file), "Raw frame this was calibrated from")
                frame.ccd.header["RAWDIR"] = (
                    os.path.dirname(os.path.abspath(sci_file)),
                    "Directory holding RAWFILE at calibration time")
                # The delivered PSF of this frame. Measured here because this is
                # where the per-frame work already happens, and recorded so that
                # phase 2 can weight and reject on it and phase 3 can size its
                # apertures from a real number instead of a config default.
                if config.phase1.steps.enabled("measure_fwhm"):
                    _record_frame_fwhm(frame.ccd.header, frame.ccd.data, config,
                                       instrument)
                frame.ccd.header["CALVERS"] = (CALVERS, "Calibration vintage")
                write_mef(
                    out_path,
                    sci=frame.ccd.data,
                    err=frame.ccd.uncertainty.array if frame.ccd.uncertainty is not None else None,
                    dq=frame.dq,
                    header=frame.ccd.header,
                    history="PHASE 1: ISR (bias/dark/flat, CR, gain); SCI/ERR/DQ in electrons",
                )

    logger.info("Phase 1 complete. Calibrated SCI/ERR/DQ frames written to output.")
