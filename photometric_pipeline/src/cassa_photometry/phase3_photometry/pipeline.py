"""Phase 3 orchestrator: batch photometry, flux calibration, and catalogs.

For each master frame it detects the science band from the header, computes a
filter-wise zero point (with uncertainty), writes a flux-calibrated image and an
error-carrying source catalog. Products go to ``<work>/phase3`` by default, or to
an explicit ``outdir``.
"""

import os

from astropy.io import fits

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import open_fits
from cassa_photometry.instruments import get_profile
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import sibling_phase_dir
from cassa_photometry.phase3_photometry.engine import UniversalPhotometryEngine
from cassa_photometry.steps import resolve_plan, resolved_names
from cassa_photometry.targets import SUPERNOVA, TargetRegistry


def _custom_action(toggles, name, engine, path, out_dir, base, logger):
    """Bind a user-supplied phase 3 step to this frame, or None if not custom.

    A custom step is called with keyword arguments so it can accept only what it
    cares about (``def defringe(*, engine, **_)``) and keep working when this
    call site later grows another. That is what stops a plugin from breaking on
    an upgrade it had no part in.
    """
    for step in resolve_plan("phase3", toggles):
        if step.name == name and step.function is not None:
            logger.info("%s: running custom step %r.", base, name)
            return lambda: step.function(
                engine=engine, path=path, outdir=out_dir, base=base,
                logger=logger,
            )
    return None


#: Returned by :func:`detect_band` for a filter that cannot be calibrated against
#: a broadband reference catalog (narrowband, or a blocked wheel position).
UNCALIBRATED = object()


def detect_band(file_path, default="R", instrument=None, config=None):
    """Return the science band (R/G/B/V/I) for a file, header first.

    Returns :data:`UNCALIBRATED` when the profile reports that this filter has no
    broadband counterpart, which is distinct from "the header did not say" -- the
    latter still falls through to the filename and then to ``default``.

    With neither ``instrument`` nor ``config`` the ``generic`` profile is used,
    which is not necessarily the profile the run was configured with; the
    pipeline always passes one.
    """
    instrument = instrument or get_profile(
        getattr(config, "instrument", None), config=config
    )
    try:
        with open_fits(file_path) as hdul:
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
        config=None, logger=None, instrument=None, targets=None):
    """Run phase 3 over a single master FITS file or a directory of them.

    ``targets`` optionally names a ``targets.yaml``; without one the frames'
    own header cards decide, and a frame that declares nothing is a field --
    which is the historical behaviour.
    """
    config = config or load_config()
    own_logger = logger is None
    logger = logger or get_logger("cassa_photometry")
    instrument = instrument or get_profile(config.instrument, config=config)
    registry = TargetRegistry(targets or config.targets_file, logger=logger)
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
    if own_logger:
        logger = get_logger("cassa_photometry", run_dir=outdir)
    logger.info(f"Phase 3 output -> {outdir}")

    logger.info(f"Phase 3: {len(files)} file(s) to process.")
    for i, path in enumerate(files, 1):
        base = os.path.basename(path).replace(".fits", "")

        band = detect_band(path, default_band, instrument)
        calibratable = band is not UNCALIBRATED
        target = _resolve_target(path, registry)
        logger.info(f"[{i}/{len(files)}] {base} | band: "
                    f"{band if calibratable else 'uncalibratable filter'} "
                    f"| target: {target.describe()}")
        _announce_routing(target, logger)

        engine = UniversalPhotometryEngine(
            fwhm_estimate=fwhm, detection_threshold=threshold, config=config, logger=logger,
        )
        _reduce_one(engine, path, base, band, calibratable, outdir, config, logger)

    logger.info("Phase 3 batch processing complete.")


def _reduce_one(engine, path, base, band, calibratable, out_dir, config, logger):
    """Run phase 3's top-level steps over one frame, in the resolved order.

    Each step is attempted SEPARATELY. Sharing one try/except meant that any
    zero-point failure -- a network blip, a blocked VizieR mirror, a field with
    no matched calibrators -- threw away the catalog, the segmentation map and
    the flux-calibrated image too, leaving an empty phase3 directory. Detection
    does not depend on calibration, and the catalog is still worth having
    without it.
    """
    def _zero_point():
        if not calibratable:
            # A narrowband or blocked filter has no broadband counterpart in the
            # reference catalogs. Cross-matching one anyway yields a zero point
            # that is numerically fine and physically meaningless, so the
            # catalog is written with instrumental magnitudes instead.
            logger.warning(
                "%s: no reference catalog covers this filter; skipping the "
                "zero point. The catalog will carry instrumental magnitudes "
                "(MAG_ISO = NaN).", base,
            )
            return
        engine.calculate_local_zero_point(path, science_band=band)

    def _flux_calibration():
        # Measuring the zero point and *writing* the calibrated image are
        # separate wants: the ZP is a number the catalog needs, while
        # _fluxcal.fits is a full-size copy of the frame that a user working
        # from catalogs alone has no use for.
        if engine.zero_point is None:
            logger.warning("%s: no zero point, so no flux-calibrated image.", base)
            return
        engine.export_flux_calibrated_image(
            path, os.path.join(out_dir, f"{base}_fluxcal.fits"))

    def _catalog():
        engine.generate_full_catalog(
            path,
            os.path.join(out_dir, f"{base}_catalog.csv"),
            os.path.join(out_dir, f"{base}_segmap.fits"),
        )

    # Nested steps (aperture_correction, psf_photometry, classification) run
    # inside these and are gated by the engine itself, so they have no entry
    # here. The registry marks them immovable for exactly that reason.
    actions = {
        "zero_point": _zero_point,
        "flux_calibration": _flux_calibration,
        "catalog": _catalog,
    }
    for name in resolved_names("phase3", config.phase3.steps):
        action = actions.get(name)
        if action is None:
            action = _custom_action(config.phase3.steps, name, engine, path,
                                    out_dir, base, logger)
            if action is None:
                continue
        try:
            action()
        except Exception as exc:
            logger.error("%s: %s failed (%s).", base, name, exc)
    logger.info("Finished %s", base)


def _resolve_target(path, registry):
    """The target a master frame describes."""

    try:
        return registry.resolve(fits.getheader(path))
    except Exception:
        return registry.resolve(None)


def _announce_routing(target, logger):
    """Say how this target will be measured, including what is not built yet.

    A supernova measured as an ordinary field source is biased low by its host's
    light, by an amount that depends where in the galaxy it sits. Saying so is
    the difference between a known limitation and a wrong number nobody queried.
    """
    if target.type == SUPERNOVA:
        logger.warning(
            "%s is declared a supernova, but host removal and forced photometry "
            "are not implemented yet (planned; see docs/master-plan.md). It will "
            "be measured as an ordinary source, which sits on its host's light "
            "and is biased LOW. Do not use these magnitudes for a light curve.",
            target.name or "This target",
        )
