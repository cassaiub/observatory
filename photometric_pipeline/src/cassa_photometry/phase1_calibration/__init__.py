"""Phase 1: Instrument Signature Removal (ISR).

Seeds the per-pixel error budget and writes multi-extension calibrated frames.
"""

from cassa_photometry.phase1_calibration.pipeline import run

__all__ = ["run"]
