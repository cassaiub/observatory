"""cassa_photometry: an end-to-end photometric reduction pipeline.

Four phases, one set of conventions:

* ``phase1_calibration`` -- Instrument Signature Removal (bias/dark/flat, cosmic
  rays, gain correction). Seeds the per-pixel error budget.
* ``phase2_integration`` -- alignment, weighted stacking, WCS astrometry.
  Propagates the ``ERR`` and ``DQ`` planes through the stack and records the
  astrometric fit residual on the WCS.
* ``phase3_photometry`` -- source detection, filter-wise zero point with
  uncertainty, flux calibration, and error-carrying catalogs. Includes
  ``verify``, an independent check against APASS/Pan-STARRS/SDSS.
* ``phase4_diagnostics`` -- per-stage PSF/QA figures, a health summary, and
  ``metrics.json``.

All science images are stored as multi-extension FITS with ``SCI``/``ERR``/``DQ``
planes (see :mod:`cassa_photometry.fits_utils`), and each phase reads and writes
its own directory under one work directory (see :mod:`cassa_photometry.paths`).
"""

__version__ = "0.2.0.dev0"

from cassa_photometry.config import PipelineConfig, load_config

__all__ = ["PipelineConfig", "load_config", "__version__"]
