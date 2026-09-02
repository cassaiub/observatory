"""Console entry points for the pipeline.

Each function is wired to a ``cassa-*`` command in ``pyproject.toml``:

* ``cassa-calibrate``  -> phase 1 (ISR)
* ``cassa-integrate``  -> phase 2 (stack + WCS)
* ``cassa-photometry`` -> phase 3 (zero point + catalogs)
* ``cassa-verify``     -> verification against reference catalogs
* ``cassa-run``        -> phases 1 -> 2 -> 3, chained
"""

import argparse

from cassa_photometry import __version__
from cassa_photometry.config import load_config
from cassa_photometry.instruments import available_profiles
from cassa_photometry.logging_utils import get_logger

_INSTRUMENT_HELP = (
    "Instrument profile to reduce with. Overrides the config file. "
    f"One of: {', '.join(available_profiles())}."
)


def _add_instrument_arg(parser):
    parser.add_argument("--instrument", default=None, help=_INSTRUMENT_HELP)


def _config_from_args(args):
    """Build the config, letting an explicit --instrument win over the YAML."""
    config = load_config(getattr(args, "config", None))
    instrument = getattr(args, "instrument", None)
    if instrument:
        config.instrument = instrument
    return config


def calibrate():
    """Phase 1: Instrument Signature Removal."""
    p = argparse.ArgumentParser(prog="cassa-calibrate", description="Phase 1: ISR calibration.")
    p.add_argument("-i", "--input", required=True, help="Directory of raw FITS files.")
    p.add_argument("-o", "--output", required=True, help="Directory for calibrated frames.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    p.add_argument("--version", action="version", version=f"cassa-photometry {__version__}")
    args = p.parse_args()

    from cassa_photometry.phase1_calibration import run
    run(args.input, args.output, config=_config_from_args(args))


def integrate():
    """Phase 2: alignment, stacking, and WCS astrometry."""
    p = argparse.ArgumentParser(prog="cassa-integrate", description="Phase 2: integration + WCS.")
    p.add_argument("data_dir", help="Directory containing calibrated FITS frames (phase1).")
    p.add_argument("-o", "--output", default=None,
                   help="Output directory (default: the phase2 directory beside the input).")
    p.add_argument("--keep-temps", action="store_true", help="Keep solve-field temporary files.")
    p.add_argument("-y", "--yes", action="store_true", help="Skip the interactive confirmation.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    args = p.parse_args()

    from cassa_photometry.phase2_integration import run
    run(args.data_dir, output_dir=args.output, keep_temps=args.keep_temps,
        config=_config_from_args(args), assume_yes=args.yes)


def photometry():
    """Phase 3: zero point, flux calibration, and catalogs."""
    p = argparse.ArgumentParser(prog="cassa-photometry", description="Phase 3: photometry + catalogs.")
    p.add_argument("input_path", help="A master FITS file or a directory of them.")
    p.add_argument("--filter", default="R", help="Fallback science band (R/G/B/V/I). Default: R.")
    p.add_argument("--fwhm", type=float, default=None, help="FWHM estimate in pixels.")
    p.add_argument("--threshold", type=float, default=None, help="Detection threshold in sigma.")
    p.add_argument("--outdir", default=None, help="Output directory (default: alongside inputs).")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    args = p.parse_args()

    from cassa_photometry.phase3_photometry import run
    run(args.input_path, default_band=args.filter, fwhm=args.fwhm, threshold=args.threshold,
        outdir=args.outdir, config=_config_from_args(args))


def verify():
    """Verify photometric calibration against reference catalogs."""
    p = argparse.ArgumentParser(prog="cassa-verify", description="Verify calibration vs catalogs.")
    p.add_argument("input_path", help="A _catalog.csv file or a directory of them.")
    p.add_argument("--filter", default="rmag", help="Fallback filter band. Default: rmag.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    args = p.parse_args()

    from cassa_photometry.phase3_photometry.verify import run as verify_run
    verify_run(args.input_path, default_filter=args.filter, config=_config_from_args(args))


def diagnose():
    """Phase 4: PSF + visual diagnostics for each pipeline stage."""
    p = argparse.ArgumentParser(prog="cassa-diagnose", description="Phase 4: pipeline diagnostics.")
    p.add_argument("run_dir", nargs="?", default=None,
                   help="A Phase 2 run directory (auto-discovers master/catalog/calibrated).")
    p.add_argument("--raw", default=None, help="Raw frame or raw directory (stage_0).")
    p.add_argument("--calibrated", default=None, help="A calibrated_*.fits frame (stage_1).")
    p.add_argument("--master", default=None, help="A Master_*.fits stack (stage_2).")
    p.add_argument("--catalog", default=None, help="A *_catalog.csv (stage_3).")
    p.add_argument("--fluxcal", default=None, help="A *_fluxcal.fits (for ZP header).")
    p.add_argument("--outdir", default=None, help="Output directory for the report.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    args = p.parse_args()

    from cassa_photometry.phase4_diagnostics import run as diag_run
    diag_run(run_dir=args.run_dir, raw=args.raw, calibrated=args.calibrated,
             master=args.master, catalog=args.catalog, fluxcal=args.fluxcal,
             outdir=args.outdir, config=_config_from_args(args))


def run_all():
    """Run phases 1 -> 2 -> 3 end to end."""
    p = argparse.ArgumentParser(prog="cassa-run", description="Run phases 1-2-3 end to end.")
    p.add_argument("-i", "--input", required=True, help="Directory of raw FITS files.")
    p.add_argument("-o", "--output", required=True,
                   help="Work directory; phases write to <output>/phase1, phase2, phase3.")
    p.add_argument("--filter", default="R", help="Fallback science band for phase 3. Default: R.")
    p.add_argument("--keep-temps", action="store_true", help="Keep solve-field temporary files.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    args = p.parse_args()

    config = _config_from_args(args)
    logger = get_logger("cassa_run")

    from cassa_photometry.phase1_calibration import run as run_p1
    from cassa_photometry.phase2_integration.pipeline import IntegrationPipeline
    from cassa_photometry.phase3_photometry import run as run_p3
    from cassa_photometry.paths import phase_dir

    p1_dir = phase_dir(args.output, 1)

    logger.info("=== Phase 1: calibration ===")
    run_p1(args.input, p1_dir, config=config, logger=logger)

    logger.info("=== Phase 2: integration ===")
    pipeline = IntegrationPipeline(p1_dir, keep_temps=args.keep_temps, config=config)
    pipeline.setup()
    pipeline.execute()

    logger.info("=== Phase 3: photometry ===")
    run_p3(pipeline.run_dir, default_band=args.filter, config=config, logger=logger)
    logger.info("End-to-end run complete.")
