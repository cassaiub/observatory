"""CASSA Observatory 8-inch instrument profile.

Imaging train
-------------
=================  ==========================================================
OTA                8-inch (200 mm) Newtonian, f/5, 1000 mm
Corrector          SharpStar 1x MPCC -- 1.0x, focal length unchanged
Camera             QHYCCD miniCAM8M -- Sony IMX585 mono, 3856 x 2180, 2.9 um
Filters            LRGB + SHO (Ha/SII/OIII) + dark cap, in the combo wheel
=================  ==========================================================

Derived: **0.598 arcsec/pixel**, field **38.4' x 21.7'**. Under 2.5" seeing that
is a ~4.2 px FWHM, so ``phase3.fwhm: 4.0`` is a better starting point for this
setup than the 3.5 default.

Focuser, mount, guide camera, guide scope, collimator, eyepiece and polarscope
are provenance, not pipeline inputs, and deliberately do not appear here.

Why the detector table is keyed by mode and gain
------------------------------------------------
The IMX585 is a CMOS sensor with a dual conversion gain (HCG/LCG). Both the
system gain in e-/ADU and the read noise are functions of the *readout mode* and
the *gain setting*, and the two modes differ sharply: at the same gain number,
LCG holds ~38 ke- of full well at ~3.5 e- of read noise while HCG holds ~9 ke- at
~1.3 e-. A single scalar -- all a CCD-era profile needs -- is simply wrong here,
so :meth:`get_gain` and :meth:`get_read_noise` interpolate the curve for this
frame's mode. The base class already passes the whole header to both, so this
needs no API change.

The mode comes from ``READOUTM`` when the header carries it. Without it the mode
is *inferred* from the gain setting against :data:`HCG_SWITCH_GAIN`, which is
right for a camera left on auto and a guess for one driven manually -- so write
``READOUTM`` at acquisition and the ambiguity disappears.

.. warning::

   ``DETECTOR_CURVES`` below is **provisional** -- published figures for the
   IMX585, not measurements of this camera. It is a placeholder so the pipeline
   degrades sensibly, not a calibration. The error budget carried through every
   phase depends on these numbers.

   Measure them with the sibling ``cassa-camchar`` package, which exists for
   exactly this and already groups its frames by gain setting::

       cassa-camchar --data /data/camchar_campaign --report

   then load the result per mode::

       Cassa8InchProfile.set_detector_curve({0: (0.61, 3.42), ...}, mode="LCG")

   Until then the honest path is to write ``EGAIN`` and ``READNOIS`` into the
   headers at acquisition: the header always wins over these curves.
"""

from bisect import bisect_left

from cassa_photometry.instruments.base import InstrumentProfile


class Cassa8InchProfile(InstrumentProfile):
    """CASSA 8-inch f/5 Newtonian with the QHY miniCAM8M (IMX585 mono)."""

    # --- Optics (documented; the header's FOCALLEN/XPIXSZ still drive phase 2) --
    APERTURE_MM = 200.0
    FOCAL_LENGTH_MM = 1000.0
    PIXEL_SIZE_UM = 2.9
    DETECTOR_SHAPE = (2180, 3856)   # (ny, nx)

    #: Readout mode -> {gain setting: (system gain [e-/ADU], read noise [e-])},
    #: for the 16-bit output.
    #:
    #: PROVISIONAL -- see the module docstring. The IMX585 has a dual conversion
    #: gain, and the two modes are genuinely different detectors: at the same
    #: gain number LCG has ~38 ke- of full well and ~3.5 e- of read noise, while
    #: HCG has ~9 ke- and ~1.3 e-. They therefore get one curve each rather than
    #: one curve with a discontinuity in it.
    DETECTOR_CURVES = {
        "LCG": {
            0:  (0.58, 3.5),
            29: (0.15, 3.2),
        },
        "HCG": {
            30:  (0.14, 1.3),
            100: (0.043, 0.9),
            200: (0.014, 0.7),
        },
    }

    #: Gain setting at or above which the camera engages HCG on its own.
    #:
    #: Only consulted when the header does not state the mode. Writing
    #: ``READOUTM`` at acquisition removes the guess entirely -- and is the only
    #: way to be right if the mode is ever driven manually.
    HCG_SWITCH_GAIN = 30

    #: Mode used when neither the header nor the gain setting can say.
    DEFAULT_READ_MODE = "LCG"

    #: Header spellings that mean each mode.
    READ_MODE_ALIASES = {
        "HCG": "HCG", "HIGH": "HCG", "HIGH GAIN": "HCG",
        "HIGH CONVERSION GAIN": "HCG", "HIGHGAIN": "HCG",
        "LCG": "LCG", "LOW": "LCG", "LOW GAIN": "LCG",
        "LOW CONVERSION GAIN": "LCG", "LOWGAIN": "LCG",
        "STANDARD": "LCG", "NORMAL": "LCG",
    }

    #: False until DETECTOR_CURVES holds real cassa-camchar measurements.
    MEASURED = False

    #: Hard clip of the 16-bit output. QHY scales the 12-bit ADC up by 16, so
    #: the ceiling is full scale rather than the generic 50000 ADU fallback;
    #: backed off to stay clear of the non-linear shoulder.
    SATURATION_ADU = 58000.0

    @property
    def name(self):
        measured = "measured" if self.MEASURED else "PROVISIONAL detector curves"
        return f"CASSA 8-inch f/5 Newtonian + QHY miniCAM8M (IMX585) [{measured}]"

    # --- Detector -------------------------------------------------------------
    @classmethod
    def set_detector_curve(cls, curve, mode=None, measured=True):
        """Replace the detector curve, e.g. with ``cassa-camchar`` PTC results.

        ``curve`` maps gain setting -> ``(gain_e_per_adu, read_noise_e)``. Pass
        ``mode`` to replace one readout mode's curve, or omit it and pass a
        ``{mode: curve}`` mapping to replace the whole table.
        """
        if mode is None:
            cls.DETECTOR_CURVES = {m: dict(c) for m, c in curve.items()}
        else:
            curves = {m: dict(c) for m, c in cls.DETECTOR_CURVES.items()}
            curves[cls.normalize_read_mode(mode) or mode] = dict(curve)
            cls.DETECTOR_CURVES = curves
        cls.MEASURED = bool(measured)

    @classmethod
    def normalize_read_mode(cls, mode):
        """Map a header's readout-mode spelling onto ``HCG``/``LCG``, or None."""
        if mode in (None, ""):
            return None
        key = str(mode).strip().upper()
        if key in cls.DETECTOR_CURVES:
            return key
        return cls.READ_MODE_ALIASES.get(key)

    def read_mode(self, header):
        """The readout mode to use for this frame.

        The header is believed first. Failing that the mode is inferred from the
        gain setting via :data:`HCG_SWITCH_GAIN`, which is right for a camera
        left on auto and a guess for one driven manually -- hence ``READOUTM``.
        """
        stated = self.normalize_read_mode(super().get_read_mode(header))
        if stated:
            return stated
        setting = self._first_float(header, ("GAIN", "EGAIN_SETTING"))
        if setting is None:
            return self.DEFAULT_READ_MODE
        return "HCG" if setting >= self.HCG_SWITCH_GAIN else "LCG"

    def _interpolate(self, header, index):
        """Interpolate one column of this frame's mode curve at its gain setting.

        Linear between the two bracketing measured points, flat outside the
        measured range -- extrapolating a conversion-gain curve past where it was
        measured invents precision that is not there.
        """
        curve = self.DETECTOR_CURVES[self.read_mode(header)]
        settings = sorted(curve)
        setting = self._first_float(header, ("GAIN", "EGAIN_SETTING"))
        if setting is None:
            setting = settings[0]

        if setting <= settings[0]:
            return curve[settings[0]][index]
        if setting >= settings[-1]:
            return curve[settings[-1]][index]

        hi = bisect_left(settings, setting)
        if settings[hi] == setting:
            return curve[setting][index]

        lo_key, hi_key = settings[hi - 1], settings[hi]
        fraction = (setting - lo_key) / (hi_key - lo_key)
        return curve[lo_key][index] + fraction * (curve[hi_key][index] - curve[lo_key][index])

    def gain_from_header(self, header):
        """What the frame claims its gain is -- ``EGAIN``/``SYSGAIN`` only.

        ``GAIN`` on this camera is the unitless *setting*, not e-/ADU, so the
        base class's ``EGAIN -> GAIN -> SYSGAIN`` cascade would misread it by
        roughly a factor of 100. ``GAIN`` keys the detector curve instead.
        """
        return self._first_float(header, ("EGAIN", "SYSGAIN"))

    def get_gain(self, header):
        """``EGAIN`` if the header carries it, else the mode+gain-keyed curve."""
        gain = self.gain_from_header(header)
        return self._interpolate(header, 0) if gain is None else gain

    def get_read_noise(self, header):
        """``READNOIS`` if the header carries it, else the mode+gain-keyed curve."""
        read_noise = self.read_noise_from_header(header)
        return self._interpolate(header, 1) if read_noise is None else read_noise

    def get_saturation(self, header):
        """Header saturation if present, else this camera's 16-bit ceiling."""
        saturation = self.saturation_from_header(header)
        return self.SATURATION_ADU if saturation is None else saturation

    # --- Optics ---------------------------------------------------------------
    def pixel_scale_arcsec(self):
        """Plate scale in arcsec/pixel from the documented optics (0.598)."""
        return self.PIXEL_SIZE_UM / self.FOCAL_LENGTH_MM * 206.265
