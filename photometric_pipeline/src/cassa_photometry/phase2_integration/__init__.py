"""Phase 2: alignment, weighted stacking, and WCS astrometry.

Propagates the ERR and DQ planes through the stack and produces WCS-solved
multi-extension master frames.
"""

from cassa_photometry.phase2_integration.pipeline import run

__all__ = ["run"]
