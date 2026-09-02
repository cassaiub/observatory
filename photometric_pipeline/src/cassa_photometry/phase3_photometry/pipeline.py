"""Phase 3 orchestrator: batch photometry, flux calibration, and catalogs.

For each master frame it detects the science band from the header, computes a
filter-wise zero point (with uncertainty), writes a flux-calibrated image and an
error-carrying source catalog. Products go to ``<work>/phase3`` by default, or to
an explicit ``outdir``.
"""

import os

from astropy.io import fits

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import sibling_phase_dir
from cassa_photometry.instruments import get_profile
from cassa_photometry.phase3_photometry.engine import UniversalPhotometryEngine

#: Returned by :func:`detect_band` for a filter that cannot be calibrated against
#: a broadband reference catalog (narrowband, or a blocked wheel position).
UNCALIBRATED = object()


def detect_band(file_path, default="R", instrument=None):
    """Return the science band (R/G/B/V/I) for a file, header first.

    Returns :data:`UNCALIBRATED` when the profile reports that this filter has no
    broadband counterpart, which is distinct from "the header did not say" -- the
    latter still falls through to the filename and then to ``default``.
    """
    instrument = instrument or get_profile()
    try:
        with fits.open(file_path) as hdul:
            band = instrument.science_band(hdul[0].header, default=default)
        if band is None:
            return UNCALIBRATED
        if band:
            return band
    except Exception:
        pass
    # Fall back to filename tokens.
    base = os.path.basename(file_path).upper()
    for b in ("R", "V", "B", "I", "G"):
        if f"_{b}_" in base or f"-{b}-" in base:
            return b
    return default


def run(input_path, default_band="R", fwhm=None, threshold=None, outdir=None,
        config=None, logger=None, instrument=None):
    """Run phase 3 over a single master FITS file or a directory of them."""
    config = config or load_config()
    logger = logger or get_logger("cassa_photometry")
    instrument = instrument or get_profile(config.instrument)
    input_path = os.path.abspath(input_path)

    if os.path.isdir(input_path):
        files = sorted(
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(".fits") and not f.endswith(("_fluxcal.fits", "_segmap.fits"))
        )
    elif os.path.isfile(input_path) and input_path.endswith(".fits"):
        files = [input_path]
    else:
        logger.error("Input must be a .fits file or a directory containing .fits files.")
        return

    if not files:
        logger.error("No master FITS files found to process.")
        return

    # Phase 3 products go to <work>/phase3, beside the phase 2 stacks they came from.
    if outdir is None:
        anchor = input_path if os.path.isdir(input_path) else os.path.dirname(input_path)
        outdir = sibling_phase_dir(anchor, 3)
    os.makedirs(outdir, exist_ok=True)
    logger.info(f"Phase 3 output -> {outdir}")

    logger.info(f"Phase 3: {len(files)} file(s) to process.")
    for i, path in enumerate(files, 1):
        base = os.path.basename(path).replace(".fits", "")
        out_dir = outdir

        band = detect_band(path, default_band, instrument)
        calibratable = band is not UNCALIBRATED
        logger.info(f"[{i}/{len(files)}] {base} | band: "
                    f"{band if calibratable else 'uncalibratable filter'}")

        engine = UniversalPhotometryEngine(
            fwhm_estimate=fwhm, detection_threshold=threshold, config=config, logger=logger,
        )
        try:
            if calibratable:
                engine.calculate_local_zero_point(path, science_band=band)
                engine.export_flux_calibrated_image(
                    path, os.path.join(out_dir, f"{base}_fluxcal.fits"))
            else:
                # A narrowband or blocked filter has no broadband counterpart in
                # APASS/Pan-STARRS/SDSS. Cross-matching one anyway yields a zero
                # point that is numerically fine and physically meaningless, so
                # the catalog is written with instrumental magnitudes instead.
                logger.warning(
                    f"{base}: no reference catalog covers this filter; skipping the "
                    f"zero point and flux calibration. The catalog will carry "
                    f"instrumental magnitudes (Absolute_Mag = NaN)."
                )
            engine.generate_full_catalog(
                path,
                os.path.join(out_dir, f"{base}_catalog.csv"),
                os.path.join(out_dir, f"{base}_segmap.fits"),
            )
            logger.info(f"Finished {base}")
        except Exception as exc:
            logger.error(f"Failed on {base}: {exc}")

    logger.info("Phase 3 batch processing complete.")
