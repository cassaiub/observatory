"""The instrument-signature-removal engine.

Applies the standard reduction steps to each science frame. Because every frame
carries a ``StdDevUncertainty`` (seeded in :mod:`~cassa_photometry.phase1_calibration.data_models`)
and every master frame carries one too, ``ccdproc`` propagates the error plane
automatically through bias/dark subtraction, flat fielding and gain correction.
A data-quality bitmask is assembled from saturation, bad-pixel and cosmic-ray
information.
"""

from dataclasses import dataclass

import astropy.units as u
import astroscrappy
import ccdproc
import numpy as np
from astropy.nddata import CCDData

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import build_dq
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.steps import resolve_plan


@dataclass
class CalibratedFrame:
    """A reduced frame: the ``CCDData`` (electrons, with uncertainty) plus its DQ."""

    ccd: CCDData
    dq: np.ndarray


def subtract_scaled_dark(ccd, master_dark, data_exposure, logger=None):
    """Subtract a master dark, rescaled from its own exposure to the frame's.

    ``ccdproc.subtract_dark`` is handed the master's *recorded* exposure rather
    than the frame's, so a 120 s master dark applied to a 60 s frame is halved
    instead of being subtracted whole.

    Scaling is only meaningful once the bias pedestal has been removed: a raw
    dark is ``bias + rate * t``, and scaling that scales the bias with it. The
    master records whether it was bias-subtracted, and we warn rather than
    silently producing a mis-scaled pedestal.
    """
    dark_exposure = master_dark.meta.get("EXPTIME")
    if dark_exposure is None or not np.isfinite(dark_exposure) or dark_exposure <= 0:
        if logger:
            logger.warning(
                "Master dark carries no usable EXPTIME; subtracting it unscaled. "
                "This is only correct if it matches the frame exposure."
            )
        dark_exposure = data_exposure

    scale = not np.isclose(float(dark_exposure), float(data_exposure), rtol=1e-3)
    if scale and not master_dark.meta.get("BIASSUB", False):
        if logger:
            logger.warning(
                f"Rescaling a master dark that was not bias-subtracted "
                f"({dark_exposure:g}s -> {data_exposure:g}s): the bias pedestal is "
                f"scaled along with the dark current. Supply bias frames to fix this."
            )

    return ccdproc.subtract_dark(
        ccd, master_dark,
        dark_exposure=float(dark_exposure) * u.s,
        data_exposure=float(data_exposure) * u.s,
        scale=scale,
    )


def crop_master_to_frame(master, frame_shape, origin, label, logger=None):
    """Align a full-frame master calibration to a subframed science frame.

    A camera windowed to an ROI is a routine setup, not an error: the master is
    cropped to the frame's region using the ``XORGSUBF``/``YORGSUBF`` origin the
    profile reports. Returns the master unchanged when the shapes already agree,
    a cropped view when the frame sits inside it, and ``None`` when the geometry
    cannot be reconciled at all -- a binning difference, say, which no amount of
    cropping fixes.

    Works for both ``CCDData`` masters (slicing carries the uncertainty plane
    with it) and plain boolean arrays such as the bad-pixel mask.
    """
    if master is None:
        return None

    data = master.data if hasattr(master, "data") else master
    if data.shape == frame_shape:
        return master

    (ny, nx), (my, mx) = frame_shape, data.shape
    x0, y0 = origin
    if 0 <= x0 and 0 <= y0 and y0 + ny <= my and x0 + nx <= mx:
        if logger:
            logger.info(
                f"Master {label} {data.shape} cropped to the science ROI "
                f"{frame_shape} at origin ({x0}, {y0})."
            )
        return master[y0:y0 + ny, x0:x0 + nx]

    if logger:
        hint = (f"the ROI origin ({x0}, {y0}) does not place a {ny}x{nx} frame "
                f"inside a {my}x{mx} master"
                if (x0 or y0) else
                "no subframe origin in the header (XORGSUBF/YORGSUBF), so this "
                "looks like a binning mismatch -- check XBINNING/YBINNING")
        logger.error(
            f"Cannot reconcile master {label} {data.shape} with science frame "
            f"{frame_shape}: {hint}. Skipping frame."
        )
    return None


def _trim_like(mask, fits_section):
    """Apply a FITS section to a boolean mask, matching ``trim_image``."""
    if mask is None:
        return None
    try:
        from astropy.nddata.utils import extract_array  # noqa: F401  (availability check)

        body = fits_section.strip()[1:-1]
        x_range, y_range = body.split(",")
        x0, x1 = (int(v) for v in x_range.split(":"))
        y0, y1 = (int(v) for v in y_range.split(":"))
        # FITS sections are 1-based and inclusive.
        return mask[y0 - 1:y1, x0 - 1:x1]
    except Exception:
        return None


def implausible_cosmic_rays(crmask, fraction_limit):
    """Fraction of the frame claimed, when a cosmic-ray mask cannot be real.

    Returns 0.0 when the mask is plausible, so the caller reads as a guard.

    Cosmic rays arrive at a rate set by the sky, not by what the telescope is
    pointed at: a few events per cm^2 per minute, each marking a handful of
    pixels. Even a large sensor in a long exposure lands well under 0.1% of its
    pixels. A mask covering percent of a frame is therefore never cosmic rays --
    it is astroscrappy's fine-structure model meeting something it cannot model,
    which in practice is a crowded field or a large resolved object whose
    stellar peaks are exactly the sharp features it hunts.

    The distinction matters because the step does not only flag: it *replaces*
    what it flags. On a globular cluster the repair deletes real starlight from
    the science frame, and the DQ flags it leaves behind then disqualify those
    same stars as zero-point calibrators in phase 3. Both failures are silent,
    and neither is recoverable from the calibrated frame afterwards.
    """
    if crmask is None or fraction_limit is None or fraction_limit >= 1.0:
        return 0.0
    fraction = float(np.count_nonzero(crmask)) / float(crmask.size)
    return fraction if fraction > fraction_limit else 0.0


class UniversalProcessor:
    """Applies ISR to science frames, propagating uncertainty and DQ."""

    def __init__(self, master_bias=None, master_dark=None, master_flats=None,
                 bpm=None, config=None, logger=None):
        self.master_bias = master_bias
        self.master_dark = master_dark
        self.master_flats = master_flats or {}
        self.bpm = bpm  # boolean bad-pixel mask (True == bad), or None
        self.config = config or load_config()
        self.logger = logger or get_logger("cassa_calibrate")
        #: Warnings already issued, so a per-frame condition is reported once.
        self._warned = set()

    def _warn_once(self, key, message, *args):
        if key in self._warned:
            return
        self._warned.add(key)
        self.logger.warning(message, *args)

    def _detect_cosmics(self, data, meta, inmask, ccd):
        """Find and repair cosmic rays. Returns ``(mask, cleaned_data)``.

        A documented extension point: subclass and override this to use a
        different detector without editing the reduction around it. See
        ``docs/CUSTOMIZING.md``.

        The ``satlevel`` handed to astroscrappy matters. Its default is 65535,
        but by this point the frame is bias-, dark- and flat-corrected, so a
        genuinely saturated stellar core sits far below that -- it is then
        detected as a cosmic ray and *replaced*, while the ERR plane still
        describes the original pixel.
        """
        cfg = self.config.phase1
        return astroscrappy.detect_cosmics(
            data,
            inmask=inmask,
            satlevel=self._saturation_level(meta, ccd),
            gain=meta["gain"], readnoise=meta["read_noise"],
            sigclip=cfg.cr_sigclip, sigfrac=cfg.cr_sigfrac, objlim=cfg.cr_objlim,
        )

    def _saturation_level(self, meta, ccd):
        """Saturation level in the frame's current units, for cosmic-ray rejection.

        The threshold is defined in raw ADU, but by the time cosmic rays are
        searched for the pedestal has been removed, so the level a saturated
        pixel now sits at is lower by that much.
        """
        level = meta.get("saturation_adu")
        if not level:
            return None
        return float(level) - float(meta.get("bias_level", 0.0) or 0.0)

    # --- ISR stages -----------------------------------------------------------
    #
    # Each stage takes the frame state and returns it. They were inline,
    # numbered `1.` to `7.` in one 150-line method, which is why phase 1 could
    # never honour a custom order: there was nothing to reorder. Splitting them
    # changed no arithmetic -- every stage below is the original code moved, not
    # rewritten -- so a default reduction is byte-identical.

    def _stage_overscan(self, state):
        """Subtract and trim the overscan. Changes the array shape.

        The trim must keep the *science* region. Trimming to the overscan
        section instead -- which is what this used to do -- throws the image
        away and keeps the strip. A profile that reports an overscan region
        without a trim region cannot be honoured, so the correction is skipped
        and said out loud rather than applied wrongly.
        """
        meta = state.meta
        if not meta.get("overscan"):
            return state
        trim = meta.get("trim")
        if not trim:
            self.logger.warning(
                "%s declares an overscan region (%s) but no trim region "
                "(TRIMSEC/DATASEC); skipping overscan subtraction.",
                meta.get("filter", "frame"), meta["overscan"],
            )
            return state
        state.ccd = ccdproc.subtract_overscan(
            state.ccd, fits_section=meta["overscan"], median=True)
        state.ccd = ccdproc.trim_image(state.ccd, fits_section=trim)
        # The saturation mask was built on the untrimmed array.
        state.sat_mask = _trim_like(state.sat_mask, trim)
        return state

    def _stage_bias(self, state):
        """Bias subtraction (propagates uncertainty in quadrature)."""
        master = state.master("bias")
        if master is not None:
            state.ccd = ccdproc.subtract_bias(state.ccd, master)
        return state

    def _stage_dark(self, state):
        """Dark subtraction, rescaled from the master's exposure to this frame's."""
        master = state.master("dark")
        if master is not None:
            state.ccd = subtract_scaled_dark(
                state.ccd, master, state.meta["exposure"], self.logger)
        return state

    def _stage_flat(self, state):
        """Flat fielding (propagates relative uncertainty)."""
        master = state.master("flat")
        if master is not None:
            state.ccd = ccdproc.flat_correct(state.ccd, master)
        else:
            self.logger.warning(
                "No master flat for filter %s; skipping flat field.",
                state.meta["filter"])
        return state

    def _stage_cosmic_rays(self, state):
        """Find and repair cosmic rays.

        Two things the defaults get wrong on real data. astroscrappy's
        `satlevel` defaults to 65535, but by this point the frame has been
        bias-, dark- and flat-corrected, so a genuinely saturated star core sits
        far below that -- it is then detected as a cosmic ray and *replaced*,
        while the ERR plane still describes the original pixel. Passing the real
        level, and masking known-saturated pixels out of the search, keeps
        stellar cores intact.
        """
        ccd, meta = state.ccd, state.meta
        saturated = state.saturated()
        inmask = state.master("bad-pixel mask")
        if saturated is not None and saturated.shape == ccd.data.shape:
            inmask = saturated if inmask is None else (inmask | saturated)

        crmask, clean_data = self._detect_cosmics(ccd.data, meta, inmask, ccd)

        # Never "repair" a saturated pixel: it is flagged, not fixed.
        if saturated is not None and saturated.shape == ccd.data.shape:
            crmask = crmask & ~saturated
            clean_data = np.where(saturated, ccd.data, clean_data)

        # A mask no cosmic-ray rate can explain is thrown away whole rather than
        # applied. Keeping the flags without the repair is not the safer half:
        # phase 3 disqualifies any calibrator star whose aperture holds a
        # flagged pixel, so a mask that sits on the stars removes the very
        # sources the zero point needs.
        claimed = implausible_cosmic_rays(crmask, self.config.phase1.cr_max_fraction)
        if claimed:
            self._warn_once(
                "cr_veto",
                "Cosmic-ray rejection flagged %.1f%% of a frame, far more than "
                "any cosmic-ray rate can produce. The mask is being discarded "
                "and the frames left uncleaned: this is what a crowded field or "
                "a large resolved object does to astroscrappy, and applying it "
                "would delete real starlight. Reject cosmic rays in the phase 2 "
                "stack instead, where they do not repeat between frames. Raise "
                "phase1.cr_max_fraction to accept the mask anyway.",
                100.0 * claimed,
            )
            ccd.header["CRVETO"] = (
                round(claimed, 4), "CR mask discarded; fraction of frame it claimed")
            crmask = np.zeros(ccd.data.shape, dtype=bool)
            clean_data = ccd.data

        ccd.data = clean_data
        state.crmask = crmask
        return state

    #: Plan step name -> the method that performs it. Steps absent here run
    #: elsewhere in the phase (linearity in data_models, the bad-pixel mask and
    #: the FWHM measurement in pipeline.run), which is why the registry marks
    #: those three immovable.
    ISR_STAGES = {
        "overscan": "_stage_overscan",
        "bias": "_stage_bias",
        "dark": "_stage_dark",
        "flat": "_stage_flat",
        "cosmic_rays": "_stage_cosmic_rays",
    }

    def process_science_frame(self, standard_ccd_list):
        """Reduce a list of :class:`StandardCCD` and return :class:`CalibratedFrame`."""
        cfg = self.config.phase1
        steps = cfg.steps
        plan = resolve_plan("phase1", steps)
        processed = []

        for std_ccd in standard_ccd_list:
            state = _FrameState(self, std_ccd)

            for step in plan:
                if step.function is not None:
                    # A user-supplied step. Keyword arguments so it can accept
                    # only what it cares about and keep working when this call
                    # site later grows another.
                    step.function(ccd=state.ccd, meta=state.meta,
                                  state=state, logger=self.logger)
                    continue
                method = self.ISR_STAGES.get(step.name)
                if method is None:
                    continue  # runs elsewhere in the phase
                if state.unusable:
                    break
                state = getattr(self, method)(state)

            if state.unusable:
                continue

            ccd = state.ccd
            # Gain correction (ADU -> electrons; scales data AND uncertainty).
            # Not a step: the frame must leave phase 1 in electrons whatever
            # else was skipped, or nothing downstream can read it.
            ccd = ccdproc.gain_correct(ccd, state.meta["gain"] * (u.electron / u.adu))
            ccd.header["BUNIT"] = "electron"

            if state.meta["fringe_needed"]:
                # Declared by the profile, not implemented by the pipeline. Say
                # so: a hook that silently does nothing is worse than no hook,
                # because the frame looks corrected and is not.
                self._warn_once(
                    "fringe",
                    "This instrument profile reports that %s frames fringe, but "
                    "no fringe correction is implemented. The fringe pattern "
                    "remains in the data.", state.meta.get("filter", "these"),
                )

            # Say what was not done, and in what order what was done ran, in the
            # file itself. A partially-reduced frame that looks finished is how
            # a wrong result gets published.
            skipped = steps.skipped()
            if skipped:
                ccd.header["CALSKIP"] = (",".join(skipped),
                                         "ISR steps deliberately skipped")
            ccd.header["CALPLAN"] = (",".join(s.name for s in plan),
                                     "ISR steps run, in order")

            # Assemble the data-quality bitmask.
            dq = build_dq(
                ccd.data.shape,
                saturated=state.saturated(),
                bad_pixel=state.master("bad-pixel mask"),
                cosmic_ray=state.crmask,
                no_data=~np.isfinite(ccd.data),
            )

            processed.append(CalibratedFrame(ccd=ccd, dq=dq))

        return processed


class _FrameState:
    """One frame in flight through the ISR stages.

    Exists so the stages can be separate callables at all: they used to share a
    dozen locals in one method body. The masters are cropped *lazily* because
    the crop depends on the frame's final geometry, which the overscan trim
    changes -- so it cannot happen before the plan runs, and must happen before
    the first stage that uses a master.
    """

    def __init__(self, processor, std_ccd):
        self._processor = processor
        self.ccd = std_ccd.ccd
        self.meta = std_ccd.meta
        self.sat_mask = std_ccd.sat_mask
        self._raw_sat_mask = std_ccd.sat_mask
        self.crmask = None
        self.unusable = False
        self._fitted = None

    def saturated(self):
        """The saturation mask, in the frame's current geometry."""
        return self._raw_sat_mask if self.sat_mask is None else self.sat_mask

    def master(self, label):
        """The named master, cropped to this frame. None when unavailable."""
        if self._fitted is None:
            self._fit_masters()
        return self._fitted.get(label)

    def _fit_masters(self):
        """Align the masters to this frame's geometry.

        A windowed camera produces science frames smaller than the full-frame
        calibrations; cropping to the ROI is the fix, and an unfixable mismatch
        is reported with its actual cause rather than a guess at binning.
        """
        processor = self._processor
        shape = self.ccd.data.shape
        origin = self.meta.get("subframe_origin", (0, 0))
        self._fitted = {}
        for label, master in (("bias", processor.master_bias),
                              ("dark", processor.master_dark),
                              ("flat", processor.master_flats.get(self.meta["filter"])),
                              ("bad-pixel mask", processor.bpm)):
            fitted = crop_master_to_frame(
                master, shape, origin, label, processor.logger)
            self._fitted[label] = fitted
            if master is not None and fitted is None:
                self.unusable = True
