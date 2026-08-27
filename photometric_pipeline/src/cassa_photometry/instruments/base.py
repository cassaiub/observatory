"""Abstract instrument profile.

Maps telescope-specific FITS metadata onto the standard model the pipeline
expects. Filter-name handling is implemented here (concrete) so that every
phase shares one convention; detector specifics (gain, read noise, amplifiers)
are left abstract for subclasses.
"""


class InstrumentProfile:
    """Base class for mapping telescope-specific FITS data to a standard model."""

    #: Raw FITS filter name (upper-cased) -> canonical stacking-group label.
    FILTER_MAP = {
        "RED": "Red", "GREEN": "Green", "BLUE": "Blue",
        "LUM": "L", "LUMINANCE": "L", "CLEAR": "L",
        "R": "R_Photo", "V": "V_Photo", "B": "B_Photo", "U": "U_Photo", "I": "I_Photo",
        "R-PHOTOMETRIC": "R_Photo", "V-PHOTOMETRIC": "V_Photo",
        "HALPHA": "Ha", "H-ALPHA": "Ha", "HA": "Ha",
        "HA-7NM": "Ha_7nm", "HA-3NM": "Ha_3nm",
        "S-II": "SII", "SII": "SII", "S2": "SII",
        "O-III": "OIII", "OIII": "OIII", "O3": "OIII",
    }

    #: Raw FITS filter name (upper-cased) -> strict science band for catalog routing.
    SCIENCE_BANDS = {
        "R": "R", "R-BAND": "R", "RMAG": "R", "R_SDSS": "R", "SDSS-R": "R",
        "V": "V",
        "B": "B",
        "I": "I", "I-BAND": "I", "IMAG": "I", "I_SDSS": "I", "SDSS-I": "I",
        "G": "G", "G-BAND": "G", "GMAG": "G", "G_SDSS": "G", "SDSS-G": "G",
    }

    @property
    def name(self):
        return "Generic Instrument"

    # --- Detector specifics (must be implemented by subclasses) ---------------
    def get_amplifiers(self, hdul):
        """Return the FITS extensions that contain science/calibration data."""
        raise NotImplementedError

    def get_image_type(self, header):
        """Standardise image type to 'bias', 'dark', 'flat', or 'science'."""
        raise NotImplementedError

    def get_exposure(self, header):
        raise NotImplementedError

    def get_gain(self, header):
        raise NotImplementedError

    def get_read_noise(self, header):
        raise NotImplementedError

    def get_overscan_region(self, header):
        """Return the overscan slice, or None if not applicable."""
        return None

    def needs_fringe_correction(self, header):
        return False

    # --- Filter handling (shared by all phases) ------------------------------
    def get_filter(self, header):
        """Return the raw filter name from the header, upper-cased."""
        return header.get("FILTER", "UNKNOWN").strip().upper()

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
        """Return the strict science band (R/G/B/V/I) for reference-catalog routing."""
        raw = self.get_filter(header)
        if not raw or raw in ("UNKNOWN", "NONE"):
            return default
        return self.SCIENCE_BANDS.get(raw, default)
