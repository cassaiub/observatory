"""cassa_photometry: an end-to-end photometric reduction pipeline.

Three phases, one set of conventions:

* ``phase1_calibration`` -- Instrument Signature Removal (bias/dark/flat, cosmic
  rays, gain correction). Seeds the per-pixel error budget.
* ``phase2_integration`` -- alignment, weighted stacking, WCS astrometry.
  Propagates the ``ERR`` and ``DQ`` planes through the stack.
* ``phase3_photometry`` -- source detection, filter-wise zero point with
  uncertainty, flux calibration, and error-carrying catalogs.

All science images are stored as multi-extension FITS with ``SCI``/``ERR``/``DQ``
planes (see :mod:`cassa_photometry.fits_utils`).
"""

__version__ = "0.1.0"

from cassa_photometry.config import PipelineConfig, load_config

__all__ = ["PipelineConfig", "load_config", "__version__"]
