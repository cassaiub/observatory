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
        Provides the saturation threshold; defaults are loaded if omitted.
    add_uncertainty : bool
        When True, attach a ``StdDevUncertainty`` computed from the detector
        gain and read noise (the start of the error budget).
    """
    if config is None:
        config = load_config()
    saturation = config.phase1.saturation_adu

    standardized_list = []
    with fits.open(filepath) as hdul:
        global_header = hdul[0].header
        amplifiers = instrument.get_amplifiers(hdul)

        for amp in amplifiers:
            header = global_header.copy()
            header.update(amp.header)

            meta = {
                "image_type": instrument.get_image_type(header),
                "exposure": instrument.get_exposure(header),
                "gain": instrument.get_gain(header),
                "read_noise": instrument.get_read_noise(header),
                "overscan": instrument.get_overscan_region(header),
                "fringe_needed": instrument.needs_fringe_correction(header),
                "filter": instrument.get_filter(header),
            }

            ccd = CCDData(amp.data, meta=header, unit=u.adu)

            # Saturation mask from the RAW ADU (the threshold is defined in ADU).
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
