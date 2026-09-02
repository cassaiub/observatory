"""Standardised CCD loading for phase 1.

Wraps raw FITS frames in Astropy ``CCDData`` objects and -- crucially -- seeds
the per-pixel error budget at the very start of the pipeline with
``ccdproc.create_deviation`` (Poisson + read noise). It also records a
saturation mask so saturated pixels can be flagged in the DQ plane later.
"""

import ccdproc
import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData

from cassa_photometry.config import load_config
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


def load_standardized_ccds(filepath, instrument, config=None, add_uncertainty=True):
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
    """
    if config is None:
        config = load_config()

    standardized_list = []
    with fits.open(filepath) as hdul:
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

            # Seed the error budget: sigma = sqrt(gain*ADU + readnoise^2) / gain.
            if add_uncertainty:
                ccd = ccdproc.create_deviation(
                    ccd,
                    gain=meta["gain"] * u.electron / u.adu,
                    readnoise=meta["read_noise"] * u.electron,
                    disregard_nan=True,
                )

            standardized_list.append(StandardCCD(ccd, meta, sat_mask=sat_mask))

    return standardized_list
