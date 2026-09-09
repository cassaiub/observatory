"""Standardised CCD loading for phase 1.

Wraps raw FITS frames in Astropy ``CCDData`` objects and -- crucially -- seeds
the per-pixel error budget at the very start of the pipeline with
``ccdproc.create_deviation`` (Poisson + read noise). It also records a
saturation mask so saturated pixels can be flagged in the DQ plane later.
"""

import astropy.units as u
import numpy as np
from astropy.nddata import CCDData, StdDevUncertainty

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import open_fits
from cassa_photometry.logging_utils import get_logger

# Missing detector constants are a property of the dataset, not of one frame, so
# the warning fires once per (profile, quantity) rather than once per file.
_WARNED = set()


def _warn_once(key, message):
    if key not in _WARNED:
        _WARNED.add(key)
        get_logger("cassa_calibrate").warning(message)


class StandardCCD:
    """A standardised wrapper: an Astropy ``CCDData`` plus unified metadata.

    Attributes
    ----------
    ccd : astropy.nddata.CCDData
        Frame data (ADU) with an attached ``StdDevUncertainty`` when seeded.
    meta : dict
        Unified metadata (exposure, gain, read noise, filter, ...).
    sat_mask : numpy.ndarray or None
        Boolean mask of pixels at/above the saturation level (from the raw ADU).
    """

    def __init__(self, ccd_data, meta, sat_mask=None):
        self.ccd = ccd_data
        self.meta = meta
        self.sat_mask = sat_mask


def _poisson_plus_read_noise(data, gain, read_noise, bias_level=0.0):
    """Per-pixel 1-sigma uncertainty in ADU.

    ``sqrt(signal_e + read_noise_e^2) / gain``, where ``signal_e`` is the
    *collected* charge -- counts above the bias pedestal, floored at zero
    because a negative measurement still carries at least the read noise.
    """
    signal_e = np.clip(
        (np.asarray(data, dtype=float) - float(bias_level or 0.0)) * float(gain), 0.0, None
    )
    return np.sqrt(signal_e + float(read_noise) ** 2) / float(gain)


def load_standardized_ccds(filepath, instrument, config=None, add_uncertainty=True,
                           bias_level=0.0):
    """Read a FITS file and return a list of :class:`StandardCCD`.

    Always returns a list (length 1 for a single-chip camera, length N for a
    multi-extension file).

    Parameters
    ----------
    filepath : str
        Path to the raw FITS file.
    instrument : InstrumentProfile
        Profile used to extract metadata (gain, read noise, filter, ...).
    config : PipelineConfig, optional
        Provides the fallback saturation threshold used when ``instrument`` does
        not know its own; defaults are loaded if omitted.
    add_uncertainty : bool
        When True, attach a ``StdDevUncertainty`` computed from the detector
        gain and read noise (the start of the error budget).
    bias_level : float
        Bias pedestal in ADU, subtracted before computing the Poisson term. The
        pedestal is an electronic offset, not collected charge, so counting it
        as signal inflates every pixel's uncertainty. Zero when no master bias
        is known yet, which is the case for the calibration frames themselves.
    """
    if config is None:
        config = load_config()

    standardized_list = []
    with open_fits(filepath) as hdul:
        global_header = hdul[0].header
        amplifiers = instrument.get_amplifiers(hdul)

        for amp in amplifiers:
            header = global_header.copy()
            header.update(amp.header)

            # Neither the header nor the profile is obliged to know the detector
            # constants. When both are silent we fall back to the configured
            # last resort rather than crashing, but say so out loud: everything
            # downstream of here treats the error budget as physical.
            gain = instrument.get_gain(header)
            if gain is None:
                gain = config.phase1.fallback_gain
                _warn_once(
                    (instrument.name, "gain"),
                    f"No gain in the header (EGAIN/GAIN/SYSGAIN) and "
                    f"'{instrument.name}' has none for this detector; assuming "
                    f"{gain} e-/ADU. Every uncertainty downstream is an estimate "
                    f"until EGAIN is written at acquisition or the profile is measured.",
                )
            read_noise = instrument.get_read_noise(header)
            if read_noise is None:
                read_noise = config.phase1.fallback_read_noise
                _warn_once(
                    (instrument.name, "read_noise"),
                    f"No read noise in the header (READNOIS/RDNOISE/E-NOISE) and "
                    f"'{instrument.name}' has none for this detector; assuming "
                    f"{read_noise} e-.",
                )

            meta = {
                "image_type": instrument.get_image_type(header),
                "exposure": instrument.get_exposure(header),
                "gain": gain,
                "read_noise": read_noise,
                "overscan": instrument.get_overscan_region(header),
                "trim": instrument.get_trim_region(header),
                "fringe_needed": instrument.needs_fringe_correction(header),
                "filter": instrument.get_filter(header),
                # ROI origin, so full-frame masters can be cropped to a
                # windowed science frame rather than rejected against it.
                "subframe_origin": instrument.get_subframe_origin(header),
            }

            ccd = CCDData(amp.data, meta=header, unit=u.adu)

            # Saturation mask from the RAW ADU (the threshold is defined in ADU).
            # The instrument gets first say -- it may know its own full well --
            # and the config value is the fallback when it does not.
            saturation = instrument.get_saturation(header)
            if saturation is None:
                saturation = config.phase1.saturation_adu
            sat_mask = ccd.data >= saturation
            # Recorded so later steps can reason about the level rather than
            # re-deriving it: cosmic-ray rejection needs it in post-bias units.
            meta["saturation_adu"] = float(saturation)
            meta["bias_level"] = float(bias_level or 0.0)

            # Non-linearity, in raw ADU and before anything else touches the
            # data. Only does something when the profile supplies a curve; the
            # default is the identity, and the frame records which it got.
            if config.phase1.steps.enabled("linearity"):
                corrected = instrument.apply_linearity(ccd.data, header)
                linearised = corrected is not ccd.data
                if linearised:
                    ccd.data = np.asarray(corrected, dtype=ccd.data.dtype)
                meta["linearised"] = bool(linearised)

            # Seed the error budget: sigma = sqrt(signal*gain + readnoise^2)/gain.
            #
            # `signal` must be the counts the detector actually collected, which
            # is the raw value minus the bias pedestal -- Poisson noise comes
            # from photons and dark current, not from an electronic offset. A
            # typical 500-2000 ADU pedestal treated as signal injects tens of
            # electrons of fictitious noise into *every* pixel, swamping the read
            # noise on a sky-limited frame.
            if add_uncertainty:
                ccd.uncertainty = StdDevUncertainty(
                    _poisson_plus_read_noise(
                        ccd.data, meta["gain"], meta["read_noise"], bias_level
                    )
                )

            standardized_list.append(StandardCCD(ccd, meta, sat_mask=sat_mask))

    return standardized_list
