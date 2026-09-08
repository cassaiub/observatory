"""Instrument profiles: the single source of truth for detector metadata.

:class:`InstrumentProfile` defines the interface every phase relies on (gain,
read noise, pixel scale, filter naming, image classification) *and* serves as
the working ``generic`` profile for any setup that writes standard FITS
keywords. A concrete profile adds what a header cannot carry -- a CMOS
conversion-gain curve, fringe susceptibility, filter-wheel conventions.

The pipeline ships only the setups the CASSA Observatory operates. Describing
another one does not require a change here: see ``docs/CUSTOMIZING.md`` for the
``detector:`` config block, the ``instrument_module:`` local-file route, and the
entry-point group for a profile shipped by another package.

Phases select one by name through :func:`get_profile`; they never construct a
profile class directly.
"""

from cassa_photometry.instruments.base import InstrumentProfile
from cassa_photometry.instruments.cassa import Cassa8InchProfile
from cassa_photometry.instruments.overrides import with_detector_overrides
from cassa_photometry.instruments.registry import (
    ENTRY_POINT_GROUP,
    PROFILES,
    available_profiles,
    get_profile,
)

__all__ = [
    "InstrumentProfile",
    "Cassa8InchProfile",
    "PROFILES",
    "ENTRY_POINT_GROUP",
    "available_profiles",
    "get_profile",
    "with_detector_overrides",
]
