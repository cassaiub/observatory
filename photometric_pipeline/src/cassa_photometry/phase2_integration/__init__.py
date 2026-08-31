"""Phase 2: alignment, weighted stacking, and WCS astrometry.

Propagates the ERR and DQ planes through the stack and produces WCS-solved
multi-extension master frames in ``<work>/phase2``. The astrometric fit residual
is recorded alongside the WCS (``CRDER1``/``CRDER2``, ``ASTRMS``, ``ASTNSTAR``),
so the solution carries an uncertainty just as the data carries an ERR plane.
"""

from cassa_photometry.phase2_integration.pipeline import run

__all__ = ["run"]
