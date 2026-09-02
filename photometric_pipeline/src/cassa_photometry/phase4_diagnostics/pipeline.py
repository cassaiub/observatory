"""Phase 4 orchestrator: run stage diagnostics and assemble a report.

Given the Phase 2 directory (and optionally a raw frame/dir) it auto-discovers
the products of each stage -- calibrated frames from the sibling ``phase1``
directory, catalogs from ``phase3`` -- runs the matching diagnostic, and writes a
multi-page ``diagnostics_report.pdf``, per-stage PNGs, a ``metrics.json`` and a
pipeline-health summary into ``<work>/phase4``. Stages whose products are absent
are skipped cleanly, and the older nested layout still resolves.
"""

import os
import glob
import json

import numpy as np

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import find_phase_dir, sibling_phase_dir
from cassa_photometry.fits_utils import read_mef
from cassa_photometry.instruments import get_profile
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
        for f in sorted(glob.glob(os.path.join(raw_path, "*.f*t*"))):
            try:
                _, _, _, header = read_mef(f)
            except Exception:
                continue
            if instrument.get_image_type(header) == "science" \
               and (not filt or instrument.get_filter(header) == filt):
                return f
    return None


def _resolve(run_dir, raw, calibrated, master, catalog, fluxcal, instrument):
    """Resolve one coherent (raw -> calibrated -> master -> catalog) set."""
    if master is None and run_dir:
        master = _first_master(run_dir)

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
    if catalog is None and base and os.path.exists(base + "_catalog.csv"):
        catalog = base + "_catalog.csv"
    if catalog is None:
        for d in (p3_dir, run_dir):
            if not d:
                continue
            cands = sorted(glob.glob(os.path.join(d, "*_catalog.csv")))
            if cands:
                catalog = cands[0]
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
    logger = logger or get_logger("cassa_diagnose")
    instrument = instrument or get_profile(config.instrument)

    r = _resolve(run_dir, raw, calibrated, master, catalog, fluxcal, instrument)
    if outdir is None:
        anchor = run_dir or (os.path.dirname(r["master"]) if r["master"] else None)
        outdir = sibling_phase_dir(anchor, 4) if anchor else "."
    diag_dir = outdir
    os.makedirs(diag_dir, exist_ok=True)
    logger.info(f"Diagnostics -> {diag_dir} (filter '{r['filter'] or '?'}')")

    report = Report(os.path.join(diag_dir, "diagnostics_report.pdf"))
    metrics, health = {}, []

    # (stage key, png stem, callable producing (metrics, figures))
    plan = [
        ("stage_0_raw", r["raw"] is not None,
         lambda: stages.diagnose_raw(r["raw"], config, logger)),
        ("stage_1_calibrated", r["calibrated"] is not None,
         lambda: stages.diagnose_calibrated(r["calibrated"], config, logger, raw_path=r["raw"])),
        ("stage_2_master", r["master"] is not None,
         lambda: stages.diagnose_master(r["master"], config, logger, single_calibrated_path=r["calibrated"])),
        ("stage_3_photometry", r["catalog"] is not None,
         lambda: stages.diagnose_photometry(r["catalog"], config, logger, fluxcal_fits=r["fluxcal"])),
    ]

    for key, present, fn in plan:
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

    report.close()

    with open(os.path.join(diag_dir, "metrics.json"), "w") as fh:
        json.dump(_jsonify({"inputs": r, "metrics": metrics}), fh, indent=2)

    logger.info("=" * 46)
    logger.info("PIPELINE HEALTH SUMMARY")
    for key, status, note in health:
        logger.info(f"  {key:<22} {status:<5} {note}")
    logger.info(f"Report: {os.path.join(diag_dir, 'diagnostics_report.pdf')}")
    logger.info("=" * 46)
    return metrics
