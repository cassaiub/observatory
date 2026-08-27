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

            # 1. Overscan subtraction
            if meta["overscan"]:
                ccd = ccdproc.subtract_overscan(ccd, overscan=ccd[:, meta["overscan"]], median=True)
                ccd = ccdproc.trim_image(ccd, fits_section=meta["overscan"])

            # 2. Bias subtraction (propagates uncertainty in quadrature)
            if self.master_bias is not None:
                if ccd.shape != self.master_bias.shape:
                    self.logger.error(
                        f"Dimension mismatch: science {ccd.shape} vs master bias "
                        f"{self.master_bias.shape}. Check binning. Skipping frame."
                    )
                    continue
                ccd = ccdproc.subtract_bias(ccd, self.master_bias)

            # 3. Dark subtraction (scaled to the science exposure)
            if self.master_dark is not None:
                ccd = ccdproc.subtract_dark(
                    ccd, self.master_dark,
                    dark_exposure=meta["exposure"] * u.s,
                    data_exposure=meta["exposure"] * u.s,
                    scale=False,
                )

            # 4. Flat fielding (propagates relative uncertainty)
            matching_flat = self.master_flats.get(meta["filter"])
            if matching_flat is not None:
                ccd = ccdproc.flat_correct(ccd, matching_flat)
            else:
                self.logger.warning(f"No master flat for filter {meta['filter']}; skipping flat field.")

            # 5. Cosmic-ray rejection
            crmask, clean_data = astroscrappy.detect_cosmics(
                ccd.data,
                inmask=(self.bpm if self.bpm is not None else None),
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
                bad_pixel=self.bpm,
                cosmic_ray=crmask,
                no_data=~np.isfinite(ccd.data),
            )

            processed.append(CalibratedFrame(ccd=ccd, dq=dq))

        return processed
