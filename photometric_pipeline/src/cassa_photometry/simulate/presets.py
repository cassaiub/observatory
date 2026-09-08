"""Dataset presets and the night model that drives them.

A preset is a complete observing run: a target, a filter set, a number of
exposures, and the calibration frames to reduce them with. Everything is derived
from one seed, so ``cassa-simulate --preset workshop`` produces bit-identical
frames on every machine -- which is why the workshop dataset can be distributed
as a command rather than a download.
"""

import numpy as np

from cassa_photometry.simulate.frames import airmass_from_altitude

#: Effective wavelength per band, nm. Used to scale seeing from the DIMM's
#: reference wavelength to the observing band.
BAND_WAVELENGTH = {"B": 445.0, "V": 551.0, "R": 658.0, "I": 806.0}

#: The CASSA 8-inch, as the simulator models it.
CASSA_8INCH = {
    "telescope": "CASSA-8IN-F5",
    "camera": "QHY miniCAM8M",
    "focal_length_mm": 1000.0,
    "pixel_size_um": 2.9,
    "pixel_scale": 2.9 / 1000.0 * 206.265,   # 0.598 arcsec/pixel
    "gain_setting": 26,
    "read_mode": "HCG",
    "ccd_temp": -10.0,
    "site_lat": 23.7806,      # Independent University, Bangladesh (Dhaka)
    "site_long": 90.4074,
    "optics_fwhm_arcsec": 1.0,   # focus + collimation + optical quality
}

PRESETS = {
    # Small and fast: what the test suite runs against. Deliberately below the
    # field size the shipped astrometry indexes can solve, so tests that need a
    # WCS supply one rather than paying for a solve.
    "tiny": {
        "shape": (768, 768),
        "field": "ngc7331",
        "bands": ["V"],
        "n_science": 3,
        "exposure": 60.0,
        "nights": 1,
        "n_bias": 3,
        "n_dark": 3,
        "n_flat": 3,
        "n_faint": 400,
        "galaxy": False,
    },
    # One night of NGC 7331 in B/V/R -- the workshop dataset.
    #
    # 2048 x 1400 is a realistic windowed ROI of the IMX585's 3856 x 2180, and
    # the size matters for more than realism: at 0.598"/px it is 20.4' x 14.0',
    # which the standard 4204/4205 index files solve comfortably. A small crop
    # would need the finer 4200-4202 indexes, which most users will not have --
    # so a preset that only solved with an unusual index set would be a trap.
    "workshop": {
        "shape": (1400, 2048),
        "field": "ngc7331",
        "bands": ["B", "V", "R"],
        "n_science": 5,
        "exposure": 60.0,
        "nights": 1,
        "n_bias": 10,
        "n_dark": 10,
        "n_flat": 10,
        "n_faint": 1600,
        "galaxy": True,
    },
    # Three nights, for epoch grouping: the frames straddle UT midnight, so
    # binning on the UT date string splits one observing night into two.
    "multinight": {
        "shape": (1024, 1024),
        "field": "ngc7331",
        "bands": ["V", "R"],
        "n_science": 4,
        "exposure": 60.0,
        "nights": 3,
        "n_bias": 5,
        "n_dark": 5,
        "n_flat": 5,
        "n_faint": 600,
        "galaxy": True,
    },
}

#: The target every preset observes.
TARGET = {
    "object": "NGC7331",
    "name": "NGC7331",
    "type": "galaxy",
    "ra": 339.266875,
    "dec": 34.415778,
}


def delivered_fwhm_arcsec(dimm_zenith_seeing, airmass, band, guide_rms, optics_fwhm):
    """The PSF a frame actually delivers, from the DIMM's zenith seeing.

    The DIMM reports seeing corrected to the zenith at 500 nm for whatever star
    it was watching. A science frame's PSF is a different quantity, and larger:

        FWHM_atm       = seeing x X^0.6 x (lambda/500nm)^-0.2
        FWHM_delivered = sqrt(FWHM_atm^2 + FWHM_guide^2 + FWHM_optics^2)

    with ``FWHM_guide ~ 2.355 x GUIDERMS``. Sizing photometric apertures from the
    DIMM number directly would under-size them by roughly 30%, so the simulator
    models the difference explicitly and the pipeline is expected to measure the
    delivered value from the stars rather than read ``SEEING``.
    """
    wavelength = BAND_WAVELENGTH.get(band, 550.0)
    atmospheric = dimm_zenith_seeing * airmass**0.6 * (wavelength / 500.0) ** -0.2
    guiding = 2.3548 * guide_rms
    return float(np.sqrt(atmospheric**2 + guiding**2 + optics_fwhm**2))


class NightModel:
    """Conditions through one night: altitude, seeing, transparency.

    Each varies the way it does in practice -- the target climbs and sets, seeing
    wanders as a correlated random walk rather than independent draws, and
    transparency drifts with occasional thin cloud. Correlated variation matters:
    independent per-frame noise averages away in a stack, while a drift does not,
    and it is the drift that a per-epoch zero point has to absorb.
    """

    def __init__(self, rng, median_seeing=1.8, night_index=0):
        self.rng = rng
        self.median_seeing = float(median_seeing)
        self.night_index = int(night_index)
        self._seeing = median_seeing
        self._transparency = 1.0

    def step(self, fraction_through_night):
        """Conditions at one point in the night, ``fraction`` in [0, 1]."""
        # The target transits mid-night: altitude peaks and falls away.
        altitude = 30.0 + 45.0 * np.sin(np.pi * fraction_through_night)
        airmass = float(airmass_from_altitude(altitude))

        # Seeing: mean-reverting random walk, never unphysically good.
        self._seeing += 0.25 * (self.median_seeing - self._seeing)
        self._seeing += self.rng.normal(0.0, 0.18)
        self._seeing = float(np.clip(self._seeing, 0.7, 5.0))

        # Transparency: slow drift plus the occasional thin cloud.
        self._transparency += 0.3 * (1.0 - self._transparency)
        self._transparency *= 1.0 + self.rng.normal(0.0, 0.02)
        if self.rng.random() < 0.08:
            self._transparency *= self.rng.uniform(0.75, 0.95)
        self._transparency = float(np.clip(self._transparency, 0.4, 1.05))

        return {
            "altitude": altitude,
            "airmass": airmass,
            "dimm_seeing": self._seeing,
            "dimm_seeing_err": float(abs(self.rng.normal(0.06, 0.02))),
            "guide_rms": float(abs(self.rng.normal(0.42, 0.08))),
            "transparency": self._transparency,
        }
