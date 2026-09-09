"""Console entry points for the pipeline.

Each function is wired to a ``cassa-*`` command in ``pyproject.toml``:

* ``cassa-calibrate``  -> phase 1 (ISR)
* ``cassa-integrate``  -> phase 2 (stack + WCS)
* ``cassa-photometry`` -> phase 3 (zero point + catalogs)
* ``cassa-verify``     -> verification against reference catalogs
* ``cassa-diagnose``   -> phase 4 (per-stage QA report)
* ``cassa-run``        -> phases 1 -> 2 -> 3, chained
* ``cassa-doctor``     -> environment diagnosis
* ``cassa``            -> umbrella dispatcher for all of the above

The umbrella exists so there is one name to remember and one place for
``--help`` to list everything; the individual commands are kept because scripts
and documentation already use them.
"""

import argparse
import os
import sys

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


def _add_plan_args(parser, phases):
    """Add ``--skip`` / ``--only`` / ``--show-plan`` for the phases a command runs.

    A step may be named bare (``--skip flat``) or qualified (``--skip
    phase1.flat``). The qualified form only matters for ``cassa-run``, where two
    phases both declare ``measure_fwhm`` and a bare name would silently mean
    both -- so there, a name that is ambiguous is an error rather than a guess.
    """
    scope = phases[0] if len(phases) == 1 else "phaseN.step for a specific phase"
    parser.add_argument(
        "--skip", action="append", default=[], metavar="STEP",
        help=f"Skip a reduction step. Repeatable. Names are {scope} steps; "
             f"see --show-plan for the list.")
    parser.add_argument(
        "--only", action="append", default=[], metavar="STEP",
        help="Run ONLY these steps, skipping every other one in the phase. "
             "Repeatable.")
    parser.add_argument(
        "--show-plan", action="store_true",
        help="Print the resolved step order (with any skips applied) and exit, "
             "without reducing anything.")
    parser.set_defaults(_plan_phases=phases)


def _resolve_step_reference(config, phases, reference, flag):
    """Map a ``--skip``/``--only`` value to the (phase, step) pairs it names."""
    if "." in reference:
        phase, _, step = reference.partition(".")
        if phase not in phases:
            raise SystemExit(
                f"error: {flag} {reference}: this command does not run {phase}. "
                f"It runs: {', '.join(phases)}.")
        if step not in getattr(config, phase).steps.step_names():
            raise SystemExit(
                f"error: {flag} {reference}: {phase} has no step {step!r}. "
                f"Run with --show-plan to list them.")
        return [(phase, step)]

    matches = [p for p in phases
               if reference in getattr(config, p).steps.step_names()]
    if not matches:
        known = sorted({s for p in phases
                        for s in getattr(config, p).steps.step_names()})
        raise SystemExit(
            f"error: {flag} {reference}: no such step. "
            f"Known here: {', '.join(known)}.")
    if len(matches) > 1:
        raise SystemExit(
            f"error: {flag} {reference} is ambiguous -- "
            f"{', '.join(f'{p}.{reference}' for p in matches)} both exist. "
            f"Qualify it with the phase.")
    return [(matches[0], reference)]


def _apply_plan_args(config, args):
    """Fold ``--skip`` / ``--only`` into the config's step toggles.

    The CLI wins over the config file, which is the usual precedence: the flag
    is the more specific, more recent statement of intent.
    """
    phases = getattr(args, "_plan_phases", None)
    if not phases:
        return config

    for reference in getattr(args, "skip", []) or []:
        for phase, step in _resolve_step_reference(config, phases, reference, "--skip"):
            toggles = getattr(config, phase).steps
            if step not in toggles.exclude:
                toggles.exclude.append(step)

    only = getattr(args, "only", []) or []
    if only:
        wanted = {}
        for reference in only:
            for phase, step in _resolve_step_reference(config, phases, reference, "--only"):
                wanted.setdefault(phase, set()).add(step)
        # --only names steps in one phase; the phases it says nothing about are
        # left alone rather than emptied, so `cassa-run --only stack` still
        # calibrates. Emptying them would make the flag a foot-gun.
        for phase, keep in wanted.items():
            toggles = getattr(config, phase).steps
            for step in toggles.step_names():
                if step not in keep and step not in toggles.exclude:
                    toggles.exclude.append(step)

    if getattr(args, "skip", None) or only:
        _revalidate(config, phases)
    return config


def _revalidate(config, phases):
    """Re-check the plans after the CLI has changed them."""
    from cassa_photometry.steps import StepPlanError, resolve_plan

    for phase in phases:
        try:
            resolve_plan(phase, getattr(config, phase).steps)
        except StepPlanError as exc:
            raise SystemExit(f"error: {exc}") from None


def _maybe_show_plan(config, args):
    """Print the resolved plan and exit, when ``--show-plan`` was given."""
    if not getattr(args, "show_plan", False):
        return
    from cassa_photometry.steps import describe_plan

    for phase in getattr(args, "_plan_phases", []):
        print(describe_plan(phase, getattr(config, phase).steps))
        print()
    raise SystemExit(0)


def _add_resource_args(parser):
    """How much of the machine the run may use.

    Setting any of these also removes the interactive CPU prompt: a stated
    budget is an answer, and being asked again would make the flag useless in
    the batch jobs and notebooks where it is most wanted.
    """
    group = parser.add_argument_group("resources")
    group.add_argument("--cores", type=int, default=None, metavar="N",
                       help="Hard ceiling on cores to use. Default: unset.")
    group.add_argument("--cpu-fraction", type=float, default=None, metavar="F",
                       help="Fraction of this machine's cores to use, 0-1 "
                            "(e.g. 0.5). Default: prompt, or 0.5 when there is "
                            "no terminal.")
    group.add_argument("--max-memory-gb", type=float, default=None, metavar="GB",
                       help="Memory the run should stay inside. Advisory: the "
                            "estimate is checked against it and a run that will "
                            "not fit says so before it starts.")


def _apply_resource_args(config, args):
    """Fold the resource flags into the config, overriding any file value."""
    if getattr(args, "cores", None) is not None:
        config.phase2.max_cores = args.cores
    if getattr(args, "cpu_fraction", None) is not None:
        config.phase2.cpu_fraction = args.cpu_fraction
    if getattr(args, "max_memory_gb", None) is not None:
        config.phase2.max_memory_gb = args.max_memory_gb
    return config


def _add_version_arg(parser):
    parser.add_argument("--version", action="version", version=f"cassa-photometry {__version__}")


def _new_parser(prog, description):
    """A parser that always carries --version."""
    parser = argparse.ArgumentParser(prog=prog, description=description)
    _add_version_arg(parser)
    return parser


def _config_from_args(args):
    """Build the config, letting an explicit --instrument win over the YAML.

    An unknown profile name is reported here, against the command the user
    typed, rather than surfacing as a ``KeyError`` from deep inside a phase.
    """
    config = load_config(getattr(args, "config", None))
    instrument = getattr(args, "instrument", None)
    if instrument:
        from cassa_photometry.instruments import get_profile

        try:
            get_profile(instrument)
        except KeyError as exc:
            raise SystemExit(f"error: {exc.args[0]}") from None
        config.instrument = instrument
    _apply_plan_args(config, args)
    _maybe_show_plan(config, args)
    return config


def calibrate():
    """Phase 1: Instrument Signature Removal."""
    p = _new_parser("cassa-calibrate", "Phase 1: ISR calibration.")
    p.add_argument("-i", "--input", required=True,
                   help="Raw FITS directory. Searched recursively, so the night "
                        "tree the acquisition software writes "
                        "(<date>/BIAS|DARK|FLAT|LIGHT/<target>/) can be given as-is.")
    p.add_argument("-o", "--output", required=True, help="Directory for calibrated frames.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    _add_plan_args(p, ["phase1"])
    args = p.parse_args()

    from cassa_photometry.phase1_calibration import run
    run(args.input, args.output, config=_config_from_args(args))


def integrate():
    """Phase 2: alignment, stacking, and WCS astrometry."""
    p = _new_parser("cassa-integrate", "Phase 2: integration + WCS.")
    p.add_argument("data_dir", help="Directory containing calibrated FITS frames (phase1).")
    p.add_argument("-o", "--output", default=None,
                   help="Output directory (default: the phase2 directory beside the input).")
    p.add_argument("--keep-temps", action="store_true", help="Keep solve-field temporary files.")
    p.add_argument("-y", "--yes", action="store_true", help="Skip the interactive confirmation.")
    p.add_argument("--raw-dir", default=None,
                   help="Raw frame directory (searched recursively), for the QA PDF's "
                        "raw panel. Only needed when the raws have moved or predate "
                        "the RAWDIR header card.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    _add_resource_args(p)
    _add_plan_args(p, ["phase2"])
    args = p.parse_args()

    config = _apply_resource_args(_config_from_args(args), args)
    if args.raw_dir:
        config.phase2.raw_dir = args.raw_dir

    from cassa_photometry.phase2_integration import run
    run(args.data_dir, output_dir=args.output, keep_temps=args.keep_temps,
        config=config, assume_yes=args.yes)


def photometry():
    """Phase 3: zero point, flux calibration, and catalogs."""
    p = _new_parser("cassa-photometry", "Phase 3: photometry + catalogs.")
    p.add_argument("input_path", help="A master FITS file or a directory of them.")
    p.add_argument("--filter", default="R", help="Fallback science band (R/G/B/V/I). Default: R.")
    p.add_argument("--fwhm", type=float, default=None, help="FWHM estimate in pixels.")
    p.add_argument("--threshold", type=float, default=None, help="Detection threshold in sigma.")
    p.add_argument("--outdir", default=None, help="Output directory (default: alongside inputs).")
    p.add_argument("--offline", action="store_true",
                   help="Never query reference catalogs; use the local cache only.")
    p.add_argument("--targets", default=None, metavar="FILE",
                   help="Optional targets.yaml. Frames carrying TARGNAME/OBJTYPE "
                        "need no file.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    _add_plan_args(p, ["phase3"])
    args = p.parse_args()

    config = _config_from_args(args)
    if args.offline:
        config.phase3.offline = True

    from cassa_photometry.phase3_photometry import run
    run(args.input_path, default_band=args.filter, fwhm=args.fwhm, threshold=args.threshold,
        outdir=args.outdir, config=config, targets=args.targets)


def verify():
    """Verify photometric calibration against reference catalogs."""
    p = _new_parser("cassa-verify", "Verify calibration vs catalogs.")
    p.add_argument("input_path", help="A _catalog.csv file or a directory of them.")
    p.add_argument("--filter", default="rmag", help="Fallback filter band. Default: rmag.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    args = p.parse_args()

    from cassa_photometry.phase3_photometry.verify import run as verify_run
    verify_run(args.input_path, default_filter=args.filter, config=_config_from_args(args))


def diagnose():
    """Phase 4: PSF + visual diagnostics for each pipeline stage."""
    p = _new_parser("cassa-diagnose", "Phase 4: pipeline diagnostics.")
    p.add_argument("run_dir", nargs="?", default=None,
                   help="A Phase 2 run directory (auto-discovers master/catalog/calibrated).")
    p.add_argument("--raw", default=None,
                   help="Raw frame, or a raw directory to search recursively (stage_0).")
    p.add_argument("--calibrated", default=None, help="A calibrated_*.fits frame (stage_1).")
    p.add_argument("--master", default=None, help="A Master_*.fits stack (stage_2).")
    p.add_argument("--catalog", default=None, help="A *_catalog.csv (stage_3).")
    p.add_argument("--fluxcal", default=None, help="A *_fluxcal.fits (for ZP header).")
    p.add_argument("--outdir", default=None, help="Output directory for the report.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    _add_plan_args(p, ["phase4"])
    args = p.parse_args()

    from cassa_photometry.phase4_diagnostics import run as diag_run
    diag_run(run_dir=args.run_dir, raw=args.raw, calibrated=args.calibrated,
             master=args.master, catalog=args.catalog, fluxcal=args.fluxcal,
             outdir=args.outdir, config=_config_from_args(args))


def run_all():
    """Run phases 1 -> 2 -> 3 end to end."""
    p = _new_parser("cassa-run", "Run phases 1-2-3 end to end.")
    p.add_argument("-i", "--input", required=True,
                   help="Raw FITS directory. Searched recursively, so the night "
                        "tree the acquisition software writes "
                        "(<date>/BIAS|DARK|FLAT|LIGHT/<target>/) can be given as-is.")
    p.add_argument("-o", "--output", required=True,
                   help="Work directory; phases write to <output>/phase1, phase2, phase3.")
    p.add_argument("--filter", default="R", help="Fallback science band for phase 3. Default: R.")
    p.add_argument("--keep-temps", action="store_true", help="Keep solve-field temporary files.")
    p.add_argument("--targets", default=None, metavar="FILE",
                   help="Optional targets.yaml for phase 3.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    p.add_argument("--from", dest="from_phase", type=int, default=1, choices=(1, 2, 3),
                   help="First phase to run. Default: 1. Use it to resume a run "
                        "whose later phases failed, without recalibrating.")
    p.add_argument("--to", dest="to_phase", type=int, default=3, choices=(1, 2, 3),
                   help="Last phase to run. Default: 3.")
    _add_instrument_arg(p)
    _add_resource_args(p)
    _add_plan_args(p, ["phase1", "phase2", "phase3"])
    args = p.parse_args()

    if args.from_phase > args.to_phase:
        raise SystemExit(
            f"error: --from {args.from_phase} is after --to {args.to_phase}.")

    config = _apply_resource_args(_config_from_args(args), args)
    logger = get_logger("cassa_run")

    from cassa_photometry.paths import phase_dir
    from cassa_photometry.phase1_calibration import run as run_p1
    from cassa_photometry.phase2_integration.pipeline import IntegrationPipeline
    from cassa_photometry.phase3_photometry import run as run_p3

    p1_dir = phase_dir(args.output, 1)
    p2_dir = phase_dir(args.output, 2)

    def _runs(phase):
        return args.from_phase <= phase <= args.to_phase

    if _runs(1):
        logger.info("=== Phase 1: calibration ===")
        run_p1(args.input, p1_dir, config=config, logger=logger)
    else:
        logger.info("=== Phase 1: skipped (--from %d) ===", args.from_phase)
        if not os.path.isdir(p1_dir):
            raise SystemExit(
                f"error: --from {args.from_phase} needs phase 1 output at "
                f"{p1_dir}, which does not exist. Run phase 1 first.")

    master_dir = p2_dir
    if _runs(2):
        logger.info("=== Phase 2: integration ===")
        pipeline = IntegrationPipeline(p1_dir, keep_temps=args.keep_temps, config=config,
                                       logger=logger, assume_yes=True)
        pipeline.setup()
        pipeline.execute()
        master_dir = pipeline.run_dir
    else:
        logger.info("=== Phase 2: skipped (--from %d) ===", args.from_phase)
        if _runs(3) and not os.path.isdir(p2_dir):
            raise SystemExit(
                f"error: phase 3 needs phase 2 output at {p2_dir}, which does "
                f"not exist. Run phase 2 first.")

    if _runs(3):
        logger.info("=== Phase 3: photometry ===")
        run_p3(master_dir, default_band=args.filter, config=config, logger=logger,
               targets=args.targets)
    logger.info("Run complete (phases %d-%d).", args.from_phase, args.to_phase)


def doctor():
    """Diagnose the environment: packages, solvers, indexes, network, permissions."""
    p = _new_parser("cassa-doctor", "Check that this machine can run the pipeline.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    p.add_argument("--no-network", action="store_true",
                   help="Skip the reachability checks (for an air-gapped machine).")
    args = p.parse_args()

    from cassa_photometry import doctor as doctor_module

    checks = doctor_module.run_checks(
        config=load_config(args.config), skip_network=args.no_network
    )
    raise SystemExit(1 if doctor_module.report(checks) else 0)


def index_fetch():
    """Download the astrometry index files a field (or a directory) needs."""
    p = _new_parser(
        "cassa-index-fetch",
        "Fetch only the astrometry index files your fields need. "
        "Run it before going off-network; the pipeline otherwise fetches them "
        "as it solves.",
    )
    p.add_argument("--ra", type=float, default=None, help="Field centre RA, degrees.")
    p.add_argument("--dec", type=float, default=None, help="Field centre Dec, degrees.")
    p.add_argument("--scale", type=float, default=None, help="Plate scale, arcsec/pixel.")
    p.add_argument("--size", type=int, default=None,
                   help="Frame size in pixels (short axis).")
    p.add_argument("--from-headers", default=None, metavar="DIR",
                   help="Read pointing and scale from every FITS frame in a "
                        "directory, and fetch the union of what they need.")
    p.add_argument("--dry-run", action="store_true",
                   help="List the files and total size without downloading.")
    p.add_argument("--wide", action="store_true",
                   help="Include the finest (largest) indexes, which the pipeline "
                        "otherwise fetches only if a solve fails.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    _add_instrument_arg(p)
    args = p.parse_args()

    if args.from_headers is None and (args.ra is None or args.dec is None):
        raise SystemExit(
            "error: give --ra/--dec (with --scale and --size), or --from-headers DIR."
        )

    from cassa_photometry.index_fetch import run as fetch_run

    raise SystemExit(fetch_run(
        ra=args.ra, dec=args.dec, scale=args.scale, size=args.size,
        from_headers=args.from_headers, dry_run=args.dry_run, wide=args.wide,
        config=_config_from_args(args),
    ))


def simulate():
    """Generate a simulated CASSA 8-inch dataset with known truth."""
    from cassa_photometry.simulate import PRESETS

    p = _new_parser("cassa-simulate", "Generate simulated data with known truth.")
    p.add_argument("--preset", default="workshop", choices=sorted(PRESETS),
                   help="Dataset to generate. Default: workshop.")
    p.add_argument("-o", "--out", default="sim", help="Output directory. Default: ./sim")
    p.add_argument("--seed", type=int, default=20260902,
                   help="RNG seed. The dataset is reproducible from it.")
    p.add_argument("-c", "--config", default=None, help="Optional YAML config file.")
    args = p.parse_args()

    from cassa_photometry.simulate import run as run_sim

    run_sim(preset=args.preset, outdir=args.out, seed=args.seed,
            config=load_config(args.config))


#: Umbrella subcommand -> (entry point, one-line description).
COMMANDS = {
    "calibrate": (calibrate, "Phase 1: instrument signature removal."),
    "integrate": (integrate, "Phase 2: alignment, stacking, WCS astrometry."),
    "photometry": (photometry, "Phase 3: zero point, flux calibration, catalogs."),
    "verify": (verify, "Check the calibration against reference catalogs."),
    "diagnose": (diagnose, "Phase 4: per-stage diagnostics report."),
    "run": (run_all, "Phases 1 -> 2 -> 3, chained."),
    "doctor": (doctor, "Check that this machine can run the pipeline."),
    "index-fetch": (index_fetch, "Fetch the astrometry indexes your fields need."),
    "simulate": (simulate, "Generate simulated data with known truth."),
}


def main(argv=None):
    """The ``cassa`` umbrella command.

    Dispatches to the same functions the individual ``cassa-*`` scripts use, by
    rewriting ``sys.argv`` so each subcommand parses its own arguments exactly
    as it does when invoked directly. That keeps one implementation and one help
    text per command rather than a second, drifting copy.
    """
    argv = list(sys.argv[1:] if argv is None else argv)

    if not argv or argv[0] in ("-h", "--help", "help"):
        width = max(len(name) for name in COMMANDS)
        print(f"cassa {__version__} -- CASSA Observatory photometric pipeline\n")
        print("usage: cassa <command> [options]\n")
        print("commands:")
        for name, (_, description) in COMMANDS.items():
            print(f"  {name:<{width}}  {description}")
        print("\nRun 'cassa <command> --help' for a command's options.")
        print("Start with 'cassa doctor' if something is not working.")
        return 0

    if argv[0] in ("-V", "--version"):
        print(f"cassa-photometry {__version__}")
        return 0

    name, rest = argv[0], argv[1:]
    if name not in COMMANDS:
        print(f"cassa: unknown command {name!r}. "
              f"Choose one of: {', '.join(COMMANDS)}.", file=sys.stderr)
        return 2

    entry_point = COMMANDS[name][0]
    sys.argv = [f"cassa-{name}", *rest]
    try:
        entry_point()
    except SystemExit as exc:  # argparse and doctor both exit this way
        return exc.code if isinstance(exc.code, int) else 1
    return 0
