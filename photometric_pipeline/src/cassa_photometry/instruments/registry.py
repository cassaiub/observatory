"""Instrument-profile registry: the one place a profile is chosen by name.

Every phase resolves its profile through :func:`get_profile` from
``config.instrument``, so a run is reduced with the profile the user asked for
rather than one hardcoded at import time. Adding a setup means adding a class and
one line here.
"""

from cassa_photometry.instruments.base import InstrumentProfile
from cassa_photometry.instruments.cassa import Cassa8InchProfile
from cassa_photometry.instruments.itelescope import ITelescopeNetworkProfile

#: Config name -> profile class.
PROFILES = {
    # Standard FITS keywords only; correct for any setup that writes them.
    "generic": InstrumentProfile,
    "itelescope": ITelescopeNetworkProfile,
    "cassa8": Cassa8InchProfile,
}

DEFAULT_PROFILE = "generic"


def get_profile(name=None):
    """Return an instrument profile instance by registry name.

    ``None`` or an empty name yields the ``generic`` profile.
    """
    key = (name or DEFAULT_PROFILE).strip().lower()
    try:
        return PROFILES[key]()
    except KeyError:
        raise KeyError(
            f"Unknown instrument profile {name!r}. "
            f"Available: {', '.join(sorted(PROFILES))}."
        ) from None


def available_profiles():
    """Sorted registry names, for CLI help text and error messages."""
    return sorted(PROFILES)
