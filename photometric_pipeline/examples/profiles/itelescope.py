"""iTelescope network instrument profile -- a retired, worked example.

This shipped with the pipeline until the instrument layer was narrowed to the
setups CASSA actually operates. It is kept here because it is a realistic
example of the *second* reason to write a profile class: a value that depends
on the header (which telescope took the frame) rather than a constant, which a
``detector:`` config block cannot express.

To use it, point a config file at it -- no installation, no edit to the
package::

    # my_config.yaml
    instrument_module: examples/profiles/itelescope.py:ITelescopeNetworkProfile

Two things in here are worth reading critically rather than copying:

* ``HARDWARE`` holds **unsourced constants** with no measurement provenance and
  no binning dependence. Compare ``Cassa8InchProfile``, which carries a
  ``MEASURED`` flag saying honestly that its curves are provisional. A profile
  should say where its numbers came from.
* ``flat_proxies`` substitutes a luminance flat for a missing red one. That is
  a real convention of the network's calibration sets, and it is also
  **scientifically lossy**: L spans ~400-700 nm and R is a subset of it, so the
  wavelength-dependent part of the pixel response and any filter-specific dust
  shadows go uncorrected. The right thing about it is that the profile
  *declares* it, so the choice is visible and per-instrument rather than
  hardcoded in a phase.

Original docstring follows.


The network's frames are ordinary single-chip FITS, so all of the header parsing
comes from :class:`~cassa_photometry.instruments.base.InstrumentProfile`. What
this profile adds is the knowledge a header cannot carry: per-telescope detector
constants for frames that arrive without ``EGAIN``/``READNOIS``, which of the
telescopes fringe, and the flat-reuse convention the network's calibration sets
rely on.
"""
from cassa_photometry.instruments.base import InstrumentProfile


class ITelescopeNetworkProfile(InstrumentProfile):
    """Profile for the iTelescope hosted network (T11, T24, T32, T68, ...)."""

    #: Telescope id -> detector constants, used only when the header is silent.
    HARDWARE = {
        "T32": {"gain": 1.0, "read_noise": 9.0, "fringes": False},
        "T24": {"gain": 1.4, "read_noise": 7.0, "fringes": False},
        "T68": {"gain": 0.77, "read_noise": 1.5, "fringes": False},
        "T11": {"gain": 1.3, "read_noise": 8.0, "fringes": True},
        "DEFAULT": {"gain": 1.0, "read_noise": 10.0, "fringes": False},
    }

    #: Filters whose fringing matters on a fringe-susceptible detector.
    FRINGE_FILTERS = ("z", "y", "ha", "sii", "luminance")

    def __init__(self):
        # Kept as an instance attribute for backwards compatibility with code
        # and tests that reach into the table directly.
        self.hardware_database = dict(self.HARDWARE)

    @property
    def name(self):
        return "iTelescope Universal Network"

    def _get_telescope_id(self, header):
        """Which iTelescope took the image, as a key into the hardware table."""
        tel_string = str(
            header.get("TELESCOP", header.get("INSTRUME", "UNKNOWN"))
        ).upper()
        for tel_id in self.hardware_database:
            if tel_id in tel_string:
                return tel_id
        return "DEFAULT"

    def _hardware(self, header):
        return self.hardware_database[self._get_telescope_id(header)]

    def get_gain(self, header):
        """Header gain if present, else the telescope's tabulated value."""
        gain = super().get_gain(header)
        return self._hardware(header)["gain"] if gain is None else gain

    def get_read_noise(self, header):
        """Header read noise if present, else the telescope's tabulated value."""
        read_noise = super().get_read_noise(header)
        return self._hardware(header)["read_noise"] if read_noise is None else read_noise

    def needs_fringe_correction(self, header):
        if not self._hardware(header)["fringes"]:
            return False
        filter_name = str(header.get("FILTER", "")).lower()
        return any(f in filter_name for f in self.FRINGE_FILTERS)

    def flat_proxies(self):
        """The network's calibration sets ship luminance flats in place of red."""
        return {"RED": "LUMINANCE"}
