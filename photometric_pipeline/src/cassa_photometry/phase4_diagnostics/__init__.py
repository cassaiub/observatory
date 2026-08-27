"""Phase 4: end-to-end diagnostics.

PSF estimation and visual QA for each pipeline stage (raw -> calibrated ->
master -> photometry), assembled into a single report with a health summary.
"""

from cassa_photometry.phase4_diagnostics.pipeline import run

__all__ = ["run"]
