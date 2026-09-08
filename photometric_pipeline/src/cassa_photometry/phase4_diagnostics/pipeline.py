"""Phase 4 orchestrator: run stage diagnostics and assemble a report.

Given the Phase 2 directory (and optionally a raw frame/dir) it auto-discovers
the products of each stage -- calibrated frames from the sibling ``phase1``
directory, catalogs from ``phase3`` -- runs the matching diagnostic, and writes a
multi-page ``diagnostics_report.pdf``, per-stage PNGs, a ``metrics.json`` and a
pipeline-health summary into ``<work>/phase4``. Stages whose products are absent
are skipped cleanly, and the older nested layout still resolves.
"""

import glob
import json
import os

import numpy as np

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import read_mef
from cassa_photometry.instruments import get_profile
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import find_phase_dir, find_raw_frames, sibling_phase_dir
from cassa_photometry.phase4_diagnostics import stages
from cassa_photometry.phase4_diagnostics.plots import Report


def _header_filter(path, instrument):
    """The frame's filter, named the way the profile names it everywhere else."""
    try:
        _, _, _, header = read_mef(path)
        filt = instrument.get_filter(header)
        # "" rather than the profile's "UNKNOWN" placeholder, so that a frame
        # with no FILTER card matches anything instead of only other unknowns.
        return "" if filt == "UNKNOWN" else filt
    except Exception:
        return ""


def _first_master(run_dir):
    cands = [f for f in glob.glob(os.path.join(run_dir, "Master_*.fits"))
             if not f.endswith(("_fluxcal.fits", "_segmap.fits", "_wcs.fits"))]
    return sorted(cands)[0] if cands else None


def _match_calibrated(cal_dir, filt, instrument):
    for f in sorted(glob.glob(os.path.join(cal_dir, "calibrated_*.fits"))):
        if not filt or _header_filter(f, instrument) == filt:
            return f
    return None


def _match_raw(raw_path, filt, instrument):
    """First science frame in ``raw_path`` matching ``filt``.

    Science-vs-calibration is the profile's call, not a hardcoded IMAGETYP list,
    so a setup whose acquisition software labels frames differently still
    resolves a stage-0 frame.
    """
    if raw_path and os.path.isfile(raw_path):
        return raw_path
    if raw_path and os.path.isdir(raw_path):
        # Walks the tree: the acquisition software files a night under
        # <date>/<TYPE>/<target>, so the science frames are rarely at the top.
        for f in find_raw_frames(raw_path):
            try:
                _, _, _, header = read_mef(f)
            except Exception:
                continue
            if instrument.get_image_type(header) == "science" \
               and (not filt or instrument.get_filter(header) == filt):
                return f
    return None


def _filter_in_name(path, filt):
    """Whether a product filename carries this filter's label."""
    stem = os.path.basename(path).upper()
    token = str(filt).upper()
    return bool(token) and (f"_{token}_" in stem or f"_{token}." in stem)


def _resolve(run_dir, raw, calibrated, master, catalog, fluxcal, instrument,
             logger=None):
    """Resolve one coherent (raw -> calibrated -> master -> catalog) set."""
    if master is None and run_dir:
        master = _first_master(run_dir)

    from cassa_photometry.logging_utils import get_logger

    logger = logger or get_logger("cassa_diagnose")
    filt = _header_filter(master, instrument) if master else ""
    base = master[:-5] if master else None  # strip .fits

    if catalog is None and base and os.path.exists(base + "_catalog.csv"):
        catalog = base + "_catalog.csv"
    # Phase 3 products sit in the sibling phase3 directory (older runs kept them
    # alongside the masters, so check both).
    p3_dir = find_phase_dir(run_dir, 3) if run_dir else None
    if catalog is None and base and p3_dir:
        cand = os.path.join(p3_dir, os.path.basename(base) + "_catalog.csv")
        catalog = cand if os.path.exists(cand) else None
    # Last resort: any catalog in the phase 3 directory -- but one whose FILTER
    # matches, not simply the alphabetically first. Taking the first meant a
    # report could pair a B-band catalog with a V-band master and describe two
    # different filters in one document without saying so.
    if catalog is None:
        for d in (p3_dir, run_dir):
            if not d:
                continue
            candidates = sorted(glob.glob(os.path.join(d, "*_catalog.csv")))
            if not candidates:
                continue
            matching = [c for c in candidates if filt and _filter_in_name(c, filt)]
            if matching:
                catalog = matching[0]
            elif len(candidates) == 1:
                catalog = candidates[0]
            else:
                # Several catalogs and no way to tell which belongs to this
                # master: no stage-3 panel is better than a mismatched one.
                logger.warning(
                    "%d catalogs found and none matches the master's filter %r; "
                    "skipping the stage-3 panel rather than pairing at random.",
                    len(candidates), filt,
                )
            break
    if fluxcal is None and base:
        for cand in ([os.path.join(p3_dir, os.path.basename(base) + "_fluxcal.fits")] if p3_dir else []) \
                    + [base + "_fluxcal.fits"]:
            if os.path.exists(cand):
                fluxcal = cand
                break

    if calibrated is None:
        # Calibrated frames live in the sibling phase1 directory; fall back to the
        # parent directory for the older nested layout.
        cal_dir = find_phase_dir(run_dir, 1) if run_dir else None
        if cal_dir is None:
            cal_dir = os.path.dirname(run_dir) if run_dir else (os.path.dirname(master) if master else None)
        if cal_dir:
            calibrated = _match_calibrated(cal_dir, filt, instrument)

    raw = _match_raw(raw, filt, instrument)
    return {"raw": raw, "calibrated": calibrated, "master": master,
            "catalog": catalog, "fluxcal": fluxcal, "filter": filt}


def _jsonify(obj):
    if isinstance(obj, dict):
        return {k: _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    return obj


def run(run_dir=None, raw=None, calibrated=None, master=None, catalog=None,
        fluxcal=None, outdir=None, config=None, logger=None, instrument=None):
    """Produce the pipeline diagnostics report. Returns the metrics dict."""
    config = config or load_config()
    own_logger = logger is None
    logger = logger or get_logger("cassa_diagnose")
    instrument = instrument or get_profile(config.instrument, config=config)

    r = _resolve(run_dir, raw, calibrated, master, catalog, fluxcal, instrument,
                 logger=logger)
    if outdir is None:
        anchor = run_dir or (os.path.dirname(r["master"]) if r["master"] else None)
        outdir = sibling_phase_dir(anchor, 4) if anchor else "."
    diag_dir = outdir
    os.makedirs(diag_dir, exist_ok=True)
    if own_logger:
        logger = get_logger("cassa_diagnose", run_dir=diag_dir)
    logger.info(f"Diagnostics -> {diag_dir} (filter '{r['filter'] or '?'}')")

    report = Report(os.path.join(diag_dir, "diagnostics_report.pdf"))
    metrics, health = {}, []

    # stage key -> (product present?, callable producing (metrics, figures)).
    available = {
        "stage_0_raw": (r["raw"] is not None,
                        lambda: stages.diagnose_raw(r["raw"], config, logger)),
        "stage_1_calibrated": (r["calibrated"] is not None,
                               lambda: stages.diagnose_calibrated(r["calibrated"], config, logger, raw_path=r["raw"])),
        "stage_2_master": (r["master"] is not None,
                           lambda: stages.diagnose_master(r["master"], config, logger, single_calibrated_path=r["calibrated"])),
        "stage_3_photometry": (r["catalog"] is not None,
                               lambda: stages.diagnose_photometry(r["catalog"], config, logger, fluxcal_fits=r["fluxcal"])),
    }

    # The order comes from the resolved plan, not from this file. The four
    # stages inspect four different products and share no state, so any order
    # is valid -- which makes phase 4 the one place a user can reorder purely
    # for how they want the report to read.
    from cassa_photometry.steps import resolved_names

    running = resolved_names("phase4", config.phase4.steps)
    plan = [(key, *available[key]) for key in running]
    for key in available:
        if key not in running:
            logger.info(f"{key}: switched off (phase4.steps).")
            health.append((key, "OFF", "switched off in config"))

    # The PDF is closed even if a stage raises inside `report.add`; otherwise
    # the file is left truncated and unopenable, which looks like a worse
    # failure than the one that actually happened.
    try:
        for key, present, fn in plan:
            # "Switched off" (recorded above, when the plan was resolved) and
            # "product not found" are different facts, and the health summary
            # distinguishes them: one is a choice, the other is a gap in the run
            # being diagnosed.
            if not present:
                logger.warning(f"{key}: product not found - skipped.")
                health.append((key, "SKIP", "product not found"))
                continue
            try:
                m, figs = fn()
                for i, fig in enumerate(figs):
                    png = os.path.join(diag_dir, f"{key}{'' if i == 0 else f'_{i}'}.png")
                    report.add(fig, png_path=png)
                metrics[key] = m
                health.append((key, "OK", ""))
                logger.info(f"{key}: OK")
            except Exception as exc:
                logger.error(f"{key}: FAILED - {exc}")
                health.append((key, "FAIL", str(exc)))
    finally:
        report.close()

    # The resolved plan of EVERY phase, not just this one. metrics.json is the
    # machine-readable record of how a reduction was produced, and "which steps
    # ran" is part of that -- without it a customised run and a default one are
    # indistinguishable to anything reading the products later.
    from cassa_photometry.steps import StepPlanError, resolved_names

    plans = {}
    for phase in ("phase1", "phase2", "phase3", "phase4"):
        toggles = getattr(config, phase).steps
        try:
            plans[phase] = {"ran": resolved_names(phase, toggles),
                            "skipped": toggles.skipped()}
        except StepPlanError as exc:  # pragma: no cover - validated at load
            plans[phase] = {"error": str(exc)}

    with open(os.path.join(diag_dir, "metrics.json"), "w") as fh:
        json.dump(_jsonify({"inputs": r, "step_plans": plans,
                            "metrics": metrics}), fh, indent=2)

    logger.info("=" * 46)
    logger.info("PIPELINE HEALTH SUMMARY")
    for key, status, note in health:
        logger.info(f"  {key:<22} {status:<5} {note}")
    logger.info(f"Report: {os.path.join(diag_dir, 'diagnostics_report.pdf')}")
    logger.info("=" * 46)
    return metrics
