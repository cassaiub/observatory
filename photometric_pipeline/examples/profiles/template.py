"""Template for a new observing setup. Copy, rename, and delete what you do not need.

READ THIS FIRST: you may not need this file.
===========================================

Three ways to describe a setup, in increasing effort. Use the first that works.

**1. A config block -- no code.** If your camera simply does not record its gain
and read noise, say so in your config file and reduce with ``generic``::

    # my_config.yaml
    instrument: generic
    detector:
      gain: 1.4                 # e-/ADU
      read_noise: 7.0           # e-
      saturation_adu: 60000
      pixel_scale_arcsec: 0.40  # unbinned

That covers most rigs. Nothing below is needed.

**2. This file.** Write a profile when the answer *depends on the header* --
a CMOS conversion-gain curve, a multi-amplifier readout, an unusual filter
naming scheme. Point your config at it, with no installation::

    # my_config.yaml
    instrument_module: my_profile.py:MyObservatoryProfile

**3. An installed package.** If you want to distribute a profile, register it
under the ``cassa_photometry.instruments`` entry-point group in your own
``pyproject.toml``::

    [project.entry-points."cassa_photometry.instruments"]
    my_scope = "my_package.profiles:MyObservatoryProfile"

The rules a profile must respect
================================

* **The frame's header is the first authority.** A profile supplies what the
  header cannot say; it does not overrule what the header does say. Every
  ``get_*`` below follows the same shape: try the header, then fall back.
* **Return ``None`` when you do not know.** Do not invent a number. ``None``
  travels up to the caller, which warns by name and uses the configured
  fallback, so an assumed error budget is never mistaken for a measured one.
* **Say where your numbers came from.** A comment naming the datasheet, the
  PTC run, or the date you measured it is the difference between a profile
  someone can trust and one they have to re-derive.

Everything is optional. The base class is concrete and already handles standard
FITS keywords, so a setup with well-written headers may need to override
nothing at all.
"""

from cassa_photometry.instruments.base import InstrumentProfile


class MyObservatoryProfile(InstrumentProfile):
    """<Telescope> + <camera>: one line on the setup this describes."""

    # ------------------------------------------------------------------
    # Optics. Used for the plate-scale fallback, which phase 2 needs to
    # hand the plate solver a scale hint and to work out which astrometry
    # index files your field needs.
    # ------------------------------------------------------------------
    FOCAL_LENGTH_MM = 1000.0
    PIXEL_SIZE_UM = 2.9          # unbinned physical pitch

    # ------------------------------------------------------------------
    # Detector constants, used ONLY when the header is silent.
    # Say where these came from -- datasheet, PTC measurement, and when.
    # ------------------------------------------------------------------
    GAIN_E_PER_ADU = 1.4         # source: <datasheet / PTC run YYYY-MM-DD>
    READ_NOISE_E = 7.0           # source: <...>
    SATURATION_ADU = 60000.0     # where this detector actually clips

    @property
    def name(self):
        return "My Observatory 10-inch + ASI533MM"

    # --- Detector ------------------------------------------------------
    def get_gain(self, header):
        """Header gain if present, else this camera's tabulated value.

        Note the base class tries ``EGAIN`` before ``GAIN``: on most CMOS
        cameras ``GAIN`` is the unitless *setting*, not e-/ADU, and reading it
        as e-/ADU is wrong by roughly a factor of 100. If your camera is like
        that, do what ``Cassa8InchProfile`` does and exclude ``GAIN`` here.
        """
        gain = super().get_gain(header)
        return self.GAIN_E_PER_ADU if gain is None else gain

    def get_read_noise(self, header):
        read_noise = super().get_read_noise(header)
        return self.READ_NOISE_E if read_noise is None else read_noise

    def get_saturation(self, header):
        saturation = super().get_saturation(header)
        return self.SATURATION_ADU if saturation is None else saturation

    def pixel_scale_arcsec(self):
        """Unbinned arcsec/pixel from the optics; binning is applied by the caller."""
        return self.PIXEL_SIZE_UM / self.FOCAL_LENGTH_MM * 206.265

    # --- Optional: detector non-linearity ------------------------------
    # def apply_linearity(self, data, header):
    #     """Correct non-linearity in raw ADU. cassa-camchar measures the curve.
    #
    #     Non-linearity is 1-5% near full well on a typical CMOS sensor, and
    #     because it grows with signal it tilts the magnitude scale rather than
    #     offsetting it. Leave this out entirely until you have measured yours.
    #     """
    #     import numpy as np
    #     coefficients = [...]              # from your PTC / linearity run
    #     return np.polyval(coefficients, data)

    # --- Optional: filter naming ---------------------------------------
    # Your wheel's names -> the pipeline's stacking labels. Merged with, not
    # replacing, the standard map in the base class.
    # FILTER_MAP = {**InstrumentProfile.FILTER_MAP, "SLOAN-R": "R_Photo"}
    #
    # Your wheel's names -> catalog science bands. Getting this wrong means
    # calibrating against a catalog that does not describe your filter, which
    # produces a valid-looking and meaningless zero point.
    # SCIENCE_BANDS = {**InstrumentProfile.SCIENCE_BANDS, "SLOAN-R": "R"}
    #
    # Filters with no broadband counterpart, so photometric calibration is
    # skipped rather than faked.
    # UNCALIBRATED_FILTERS = InstrumentProfile.UNCALIBRATED_FILTERS | {"OIII-3NM"}

    # --- Optional: calibration conventions -----------------------------
    # def flat_proxies(self):
    #     """filter -> stand-in filter, when you routinely reuse a flat.
    #
    #     Think before adding one: a stand-in flat leaves the wavelength-dependent
    #     part of the pixel response uncorrected. See examples/profiles/itelescope.py.
    #     """
    #     return {"RED": "LUMINANCE"}
