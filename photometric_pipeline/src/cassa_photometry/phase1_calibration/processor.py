"""The instrument-signature-removal engine.

Applies the standard reduction steps to each science frame. Because every frame
carries a ``StdDevUncertainty`` (seeded in :mod:`~cassa_photometry.phase1_calibration.data_models`)
and every master frame carries one too, ``ccdproc`` propagates the error plane
automatically through bias/dark subtraction, flat fielding and gain correction.
A data-quality bitmask is assembled from saturation, bad-pixel and cosmic-ray
information.
"""

from dataclasses import dataclass

import numpy as np
import astropy.units as u
import ccdproc
import astroscrappy
from astropy.nddata import CCDData

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import build_dq
from cassa_photometry.logging_utils import get_logger


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

    def process_science_frame(self, standard_ccd_list):
        """Reduce a list of :class:`StandardCCD` and return :class:`CalibratedFrame`."""
        cfg = self.config.phase1
        processed = []

        for std_ccd in standard_ccd_list:
            ccd = std_ccd.ccd
            meta = std_ccd.meta

            # 1. Overscan subtraction (changes the shape, so it comes first)
            if meta["overscan"]:
                ccd = ccdproc.subtract_overscan(ccd, overscan=ccd[:, meta["overscan"]], median=True)
                ccd = ccdproc.trim_image(ccd, fits_section=meta["overscan"])

            # 1b. Align the masters to this frame's geometry. A windowed camera
            # produces science frames smaller than the full-frame calibrations;
            # cropping to the ROI is the fix, and an unfixable mismatch is
            # reported with its actual cause rather than a guess at binning.
            frame_shape = ccd.data.shape
            origin = meta.get("subframe_origin", (0, 0))
            fitted, unusable = {}, False
            for label, master in (("bias", self.master_bias),
                                  ("dark", self.master_dark),
                                  ("flat", self.master_flats.get(meta["filter"])),
                                  ("bad-pixel mask", self.bpm)):
                fitted[label] = crop_master_to_frame(
                    master, frame_shape, origin, label, self.logger)
                if master is not None and fitted[label] is None:
                    unusable = True
            if unusable:
                continue
            master_bias, master_dark = fitted["bias"], fitted["dark"]
            matching_flat, bpm = fitted["flat"], fitted["bad-pixel mask"]

            # 2. Bias subtraction (propagates uncertainty in quadrature)
            if master_bias is not None:
                ccd = ccdproc.subtract_bias(ccd, master_bias)

            # 3. Dark subtraction (rescaled from the master's exposure to this frame's)
            if master_dark is not None:
                ccd = subtract_scaled_dark(ccd, master_dark, meta["exposure"], self.logger)

            # 4. Flat fielding (propagates relative uncertainty)
            if matching_flat is not None:
                ccd = ccdproc.flat_correct(ccd, matching_flat)
            else:
                self.logger.warning(f"No master flat for filter {meta['filter']}; skipping flat field.")

            # 5. Cosmic-ray rejection
            crmask, clean_data = astroscrappy.detect_cosmics(
                ccd.data,
                inmask=(bpm if bpm is not None else None),
                gain=meta["gain"], readnoise=meta["read_noise"],
                sigclip=cfg.cr_sigclip, sigfrac=cfg.cr_sigfrac, objlim=cfg.cr_objlim,
            )
            ccd.data = clean_data

            # 6. Gain correction (ADU -> electrons; scales data AND uncertainty)
            ccd = ccdproc.gain_correct(ccd, meta["gain"] * (u.electron / u.adu))
            ccd.header["BUNIT"] = "electron"

            if meta["fringe_needed"]:
                # Fringe correction is not yet implemented for this network.
                pass

            # 7. Assemble the data-quality bitmask
            dq = build_dq(
                ccd.data.shape,
                saturated=std_ccd.sat_mask,
                bad_pixel=bpm,
                cosmic_ray=crmask,
                no_data=~np.isfinite(ccd.data),
            )

            processed.append(CalibratedFrame(ccd=ccd, dq=dq))

        return processed
