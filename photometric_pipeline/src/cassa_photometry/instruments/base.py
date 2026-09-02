"""The instrument profile: how a telescope's FITS metadata maps to the pipeline.

This class is both the interface every phase relies on *and* a working profile in
its own right. Everything here is implemented against standard FITS keywords, so
any setup whose acquisition software writes them correctly -- an amateur rig or a
professional observatory -- can be reduced with no subclass at all. It is
registered as the ``generic`` profile.

A subclass exists to supply what the *header cannot say*: detector constants for
cameras that do not record their own gain and read noise, multi-amplifier layouts,
fringe susceptibility, filter-wheel quirks. It should override the narrowest thing
that differs and inherit the rest -- see :mod:`~cassa_photometry.instruments.cassa`
for a worked example.

Detector values follow one order of authority throughout:

1. the frame's own header -- always wins;
2. the profile's hardware knowledge -- when the header is silent;
3. ``None`` -- meaning "nobody knows", which the caller reports rather than
   papering over with a fabricated number.
"""


class InstrumentProfile:
    """Maps telescope-specific FITS metadata onto the pipeline's standard model."""

    #: Raw FITS filter name (upper-cased) -> canonical stacking-group label.
    FILTER_MAP = {
        "RED": "Red", "GREEN": "Green", "BLUE": "Blue",
        "LUM": "L", "LUMINANCE": "L", "CLEAR": "L", "L": "L",
        "R": "R_Photo", "V": "V_Photo", "B": "B_Photo", "U": "U_Photo", "I": "I_Photo",
        "R-PHOTOMETRIC": "R_Photo", "V-PHOTOMETRIC": "V_Photo",
        "B-PHOTOMETRIC": "B_Photo", "I-PHOTOMETRIC": "I_Photo",
        "HALPHA": "Ha", "H-ALPHA": "Ha", "HA": "Ha",
        "HA-7NM": "Ha_7nm", "HA-3NM": "Ha_3nm",
        "S-II": "SII", "SII": "SII", "S2": "SII",
        "O-III": "OIII", "OIII": "OIII", "O3": "OIII",
    }

    #: Raw FITS filter name (upper-cased) -> strict science band for catalog routing.
    #:
    #: A filter absent from this map falls back to the caller's default band, so
    #: every filter a profile expects to see must appear here. The long
    #: ``*-PHOTOMETRIC`` forms are listed alongside the bare letters because the
    #: two naming conventions reach the pipeline from different filter wheels and
    #: an unlisted one is silently calibrated against the wrong catalog.
    SCIENCE_BANDS = {
        "R": "R", "R-BAND": "R", "RMAG": "R", "R_SDSS": "R", "SDSS-R": "R",
        "R-PHOTOMETRIC": "R", "RED": "R",
        "V": "V", "V-PHOTOMETRIC": "V",
        "B": "B", "B-BAND": "B", "B-PHOTOMETRIC": "B", "BLUE": "B",
        "I": "I", "I-BAND": "I", "IMAG": "I", "I_SDSS": "I", "SDSS-I": "I",
        "I-PHOTOMETRIC": "I",
        "G": "G", "G-BAND": "G", "GMAG": "G", "G_SDSS": "G", "SDSS-G": "G",
        "GREEN": "G",
        # A luminance/clear filter is broad (~400-700 nm); V is its closest
        # broadband counterpart. Change this one line to treat luminance as
        # unfiltered-clear instead, or drop it into UNCALIBRATED_FILTERS to
        # exclude it from photometric calibration entirely.
        "L": "V", "LUM": "V", "LUMINANCE": "V", "CLEAR": "V",
    }

    #: Filters with no meaningful broadband catalog counterpart.
    #:
    #: Narrowband frames and blocked filter-wheel positions must never be
    #: zero-pointed against a broadband reference catalog: the cross-match
    #: succeeds and returns a number, but that number is meaningless. The
    #: photometry phase skips flux calibration for these rather than emitting one.
    UNCALIBRATED_FILTERS = {
        "HA", "HALPHA", "H-ALPHA", "HA-7NM", "HA-3NM",
        "SII", "S-II", "S2", "OIII", "O-III", "O3",
        "DARK", "DARKCAP", "DARK-CAP", "BLANK", "CLOSED",
    }

    #: Header keywords consulted for each detector constant, in priority order.
    GAIN_KEYS = ("EGAIN", "GAIN", "SYSGAIN")
    READ_NOISE_KEYS = ("READNOIS", "RDNOISE", "E-NOISE")
    SATURATION_KEYS = ("SATURATE", "SATLEVEL", "FULLWELL")
    READ_MODE_KEYS = ("READOUTM", "READMODE", "READOUT")
    SUBFRAME_ORIGIN_KEYS = (("XORGSUBF", "YORGSUBF"), ("SUBFRAMX", "SUBFRAMY"))
    CALIBRATION_STATUS_KEYS = ("CALSTAT",)

    #: ``BUNIT`` values meaning the frame has already been gain-corrected.
    CALIBRATED_UNITS = {"ELECTRON", "ELECTRONS", "E-", "E", "ELECTRON/S", "E-/S"}

    @property
    def name(self):
        return "Generic Instrument"

    # --- Header helpers ------------------------------------------------------
    @staticmethod
    def _first_float(header, keys):
        """First key that parses as a float, or None if none of them do.

        Cameras write text into numeric slots more often than one would like
        (``'Mode0'`` in a read-noise card, for instance), so an unparseable value
        is skipped rather than allowed to raise.
        """
        for key in keys:
            value = header.get(key)
            if value is None:
                continue
            try:
                return float(value)
            except (TypeError, ValueError):
                continue
        return None

    # --- Detector -------------------------------------------------------------
    def get_amplifiers(self, hdul):
        """Return the FITS extensions holding science/calibration data.

        Single-chip, single-amplifier cameras -- nearly all of them -- keep their
        data in the primary HDU. Mosaic and multi-amplifier detectors override.
        """
        return [hdul[0]]

    def get_image_type(self, header):
        """Standardise ``IMAGETYP`` to 'bias', 'dark', 'flat', or 'science'.

        Matched as a lower-cased substring so the many spellings in the wild
        (``BIAS``, ``Bias Frame``, ``Flat Field``, ``Light Frame``) all land
        correctly. Anything matching none of them is science.
        """
        img_type = str(header.get("IMAGETYP", "")).lower()
        if "bias" in img_type:
            return "bias"
        if "dark" in img_type:
            return "dark"
        if "flat" in img_type:
            return "flat"
        return "science"

    def get_exposure(self, header):
        """Exposure time in seconds."""
        exposure = self._first_float(header, ("EXPTIME", "EXPOSURE"))
        return 0.0 if exposure is None else exposure

    def get_gain(self, header):
        """System gain in electrons per ADU, or None if the header is silent.

        ``EGAIN`` is tried before ``GAIN`` deliberately: on most CMOS cameras
        ``GAIN`` holds the unitless *gain setting* while ``EGAIN`` holds the
        system gain in e-/ADU. Reading the setting as e-/ADU is off by orders of
        magnitude, and silently so.
        """
        return self._first_float(header, self.GAIN_KEYS)

    def get_read_noise(self, header):
        """Read noise in electrons, or None if the header is silent."""
        return self._first_float(header, self.READ_NOISE_KEYS)

    def get_saturation(self, header):
        """Saturation level in raw ADU, or None if unknown.

        Full well is a property of the detector, not of the reduction, so a
        profile for a known camera should override this (or serve it from a
        hardware table) rather than making every run carry a config file just to
        flag saturated pixels correctly.
        """
        return self._first_float(header, self.SATURATION_KEYS)

    def get_overscan_region(self, header):
        """Return the overscan slice, or None if not applicable."""
        return None

    def needs_fringe_correction(self, header):
        return False

    def get_read_mode(self, header):
        """Detector readout mode as an upper-case string, or None.

        CMOS sensors with a dual conversion gain (the IMX585's HCG/LCG, for
        instance) have genuinely different gain and read noise at the *same* gain
        setting, so the mode is part of the detector's state, not a note.
        """
        for key in self.READ_MODE_KEYS:
            value = header.get(key)
            if value not in (None, ""):
                return str(value).strip().upper()
        return None

    def get_subframe_origin(self, header):
        """Return the ``(x0, y0)`` origin of a subframe/ROI, in image pixels.

        ``(0, 0)`` for a full frame, which is also the answer when the camera
        does not record an origin. Used to crop full-frame master calibrations
        down to a subframed science frame instead of rejecting it.
        """
        for x_key, y_key in self.SUBFRAME_ORIGIN_KEYS:
            x0 = self._first_float(header, (x_key,))
            y0 = self._first_float(header, (y_key,))
            if x0 is not None and y0 is not None:
                return int(x0), int(y0)
        return 0, 0

    def get_calibration_status(self, header):
        """``CALSTAT`` as an upper-case string; ``""`` means the frame is raw.

        The convention is one letter per step already applied (``B`` bias,
        ``D`` dark, ``F`` flat). Values meaning "nothing applied" normalise to
        ``""`` so callers only have to test for truthiness.
        """
        for key in self.CALIBRATION_STATUS_KEYS:
            value = header.get(key)
            if value in (None, ""):
                continue
            status = str(value).strip().upper()
            if status in ("", "NONE", "RAW", "N/A", "0"):
                return ""
            return status
        return ""

    def already_calibrated(self, header):
        """Why this frame looks reduced already, or None if it looks raw.

        Two independent signals, because either alone can be absent: what the
        acquisition software declares in ``CALSTAT``, and the physical unit --
        a raw frame is in ADU, and only a gain-corrected one is in electrons.
        Reducing an already-reduced frame is silently wrong rather than loud,
        so it is worth one cheap check per file.
        """
        status = self.get_calibration_status(header)
        if status:
            return f"CALSTAT={status!r}"
        unit = str(header.get("BUNIT", "")).strip().upper()
        if unit in self.CALIBRATED_UNITS:
            return f"BUNIT={unit!r} (frame is in electrons, not raw ADU)"
        return None

    # --- Calibration conventions ---------------------------------------------
    def flat_proxies(self):
        """Map of ``filter -> stand-in filter`` for missing master flats.

        Some setups routinely acquire flats in one filter and reuse them for
        another whose response is close enough. That is a property of the
        instrument, not of the reduction, so phase 1 applies whatever the profile
        declares here instead of hardcoding any one telescope's habit.
        """
        return {}

    # --- Filter handling (shared by all phases) ------------------------------
    def get_filter(self, header):
        """Return the raw filter name from the header, upper-cased."""
        return str(header.get("FILTER", "UNKNOWN")).strip().upper()

    def standardize_filter(self, header):
        """Return the canonical stacking-group filter label.

        Handles narrowband bandwidth suffixes (e.g. ``BANDWID`` -> ``_7nm``) the
        same way the integration phase groups frames.
        """
        raw = self.get_filter(header)
        bandwidth = header.get("BANDWID", "")
        suffix = f"_{bandwidth}nm" if bandwidth else ""
        if "NM" in raw and raw not in self.FILTER_MAP:
            return raw.replace("-", "_")
        return self.FILTER_MAP.get(raw, raw) + suffix

    def science_band(self, header, default="R"):
        """Return the science band (R/G/B/V/I) for reference-catalog routing.

        Returns ``None`` for a filter that has no broadband counterpart --
        narrowband, or a blocked filter-wheel position -- which tells the
        photometry phase to skip flux calibration rather than calibrate against
        a catalog that does not describe this filter.
        """
        raw = self.get_filter(header)
        if not raw or raw in ("UNKNOWN", "NONE"):
            return default
        if raw in self.UNCALIBRATED_FILTERS:
            return None
        return self.SCIENCE_BANDS.get(raw, default)
