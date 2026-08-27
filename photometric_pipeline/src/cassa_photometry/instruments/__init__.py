"""Instrument profiles: the single source of truth for detector metadata.

The base class defines the interface every phase relies on (gain, read noise,
filter naming, image classification); :class:`ITelescopeNetworkProfile` is the
concrete profile for the iTelescope network.
"""

from cassa_photometry.instruments.base import InstrumentProfile
from cassa_photometry.instruments.itelescope import ITelescopeNetworkProfile

__all__ = ["InstrumentProfile", "ITelescopeNetworkProfile"]
