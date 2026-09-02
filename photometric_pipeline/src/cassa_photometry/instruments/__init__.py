"""Instrument profiles: the single source of truth for detector metadata.

:class:`InstrumentProfile` defines the interface every phase relies on (gain,
read noise, filter naming, image classification) *and* serves as the working
``generic`` profile for any setup that writes standard FITS keywords. The
concrete profiles add what a header cannot carry -- per-telescope detector
constants, fringe susceptibility, filter-wheel conventions.

Phases select one by name through :func:`get_profile`; they never construct a
profile class directly.
"""

from cassa_photometry.instruments.base import InstrumentProfile
from cassa_photometry.instruments.cassa import Cassa8InchProfile
from cassa_photometry.instruments.itelescope import ITelescopeNetworkProfile
from cassa_photometry.instruments.registry import (
    PROFILES, available_profiles, get_profile,
)

__all__ = [
    "InstrumentProfile",
    "ITelescopeNetworkProfile",
    "Cassa8InchProfile",
    "PROFILES",
    "available_profiles",
    "get_profile",
]
