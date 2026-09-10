"""Shared setup for every workshop notebook.

Each notebook used to open with the same thirty lines of path juggling, copied
six times. One copy means a path is changed once, and that the notebooks cannot
drift apart -- which they had.

Usage, as the first cell of a notebook::

    from workshop_config import *
    show_config()
"""

import os
import sys

__all__ = [
    "WORKSHOP_DIR", "RAW_DIR", "WORK_DIR", "TRUTH_SOURCES", "TRUTH_FRAMES",
    "PHASE1_DIR", "PHASE2_DIR", "PHASE3_DIR", "PHASE4_DIR", "INSTRUMENT",
    "ASTROMETRY_INDEX_DIR",
    "show_config", "require_dataset", "raw_frames", "raw_layout",
    "load_truth", "score_against_truth",
]

# --- Paths --------------------------------------------------------------------
# Resolved relative to this file, so a notebook works whatever directory Jupyter
# was launched from. Override any of them with an environment variable.

WORKSHOP_DIR = os.path.abspath(
    os.environ.get("CASSA_WORKSHOP_DIR", os.path.join(os.path.dirname(__file__), ".."))
)
RAW_DIR = os.environ.get("RAW_DIR", os.path.join(WORKSHOP_DIR, "raw"))
WORK_DIR = os.environ.get("WORK_DIR", os.path.join(WORKSHOP_DIR, "work"))

PHASE1_DIR = os.path.join(WORK_DIR, "phase1")
PHASE2_DIR = os.path.join(WORK_DIR, "phase2")
PHASE3_DIR = os.path.join(WORK_DIR, "phase3")
PHASE4_DIR = os.path.join(WORK_DIR, "phase4")

TRUTH_SOURCES = os.path.join(WORKSHOP_DIR, "truth_sources.csv")
TRUTH_FRAMES = os.path.join(WORKSHOP_DIR, "truth_frames.csv")

#: Which instrument profile to reduce with.
#:
#: ``generic`` is correct for the shipped night: it is iTelescope data from a
#: CDK700 with an Andor DU934P CCD, and the ``cassa8`` profile describes a
#: completely different detector -- an IMX585 CMOS, whose gain curve would be
#: applied to a CCD that has nothing to do with it. Those frames carry no
#: ``EGAIN`` or ``READNOIS``, so phase 1 falls back to the configured constants
#: and warns per frame, which is the honest outcome.
#:
#: Set ``CASSA_INSTRUMENT=cassa8`` when reducing frames that really did come
#: from the CASSA 8-inch -- ``cassa-simulate`` writes those.
INSTRUMENT = os.environ.get("CASSA_INSTRUMENT", "generic")

#: A full local Astrometry.net index set, if one exists beside the checkout.
#:
#: Nothing needs this: whichever backend solves, the sky data a field needs is
#: worked out from the frame and fetched on demand into
#: ``~/.cache/cassa-photometry/``. But the pipeline's own default is
#: ``./astrometry_data`` *relative to the working directory*, which from
#: ``workshop/notebooks`` resolves to a path that does not exist -- a set that
#: ships with the checkout sits beside ``workshop/``. So point the pipeline at
#: it when it is there, unless the environment already names one, which wins: on
#: a shared machine the indexes are usually a read-only directory somebody else
#: maintains.
ASTROMETRY_INDEX_DIR = os.environ.get(
    "CASSA_ASTROMETRY_INDEX",
    os.path.join(os.path.dirname(WORKSHOP_DIR), "astrometry_data"),
)
if os.path.isdir(ASTROMETRY_INDEX_DIR):
    os.environ["CASSA_ASTROMETRY_INDEX"] = ASTROMETRY_INDEX_DIR


def _solver_summary():
    """Which plate-solving backend a run would use, and which are installed.

    Three backends can solve, and any one is enough -- the products are the
    same whichever ran, because the pipeline measures the astrometric residual
    itself rather than taking it from the solver. So the useful answer is the
    one that would actually be used, not whether a particular binary is on PATH.
    """
    try:
        from cassa_photometry.config import load_config
        from cassa_photometry.phase2_integration.solvers import available_backends

        usable = available_backends(load_config())
    except Exception as exc:  # pragma: no cover - diagnostic path only
        return f"could not be determined ({exc})"
    if not usable:
        return "NONE INSTALLED -- phase 2 cannot solve a WCS; run cassa-doctor"
    return f"{usable[0]}  (also available: {', '.join(usable[1:]) or 'none'})"


def show_config():
    """Print the resolved paths and the environment, and check the essentials."""
    import cassa_photometry

    print(f"cassa-photometry {cassa_photometry.__version__}")
    print(f"python           {sys.version.split()[0]}")
    print(f"plate solver     {_solver_summary()}")
    print()
    for label, path in (
        ("workshop", WORKSHOP_DIR), ("raw", RAW_DIR), ("work", WORK_DIR),
        ("truth", TRUTH_SOURCES), ("indexes", ASTROMETRY_INDEX_DIR),
    ):
        mark = "ok " if os.path.exists(path) else "-- "
        print(f"  [{mark}] {label:9s} {path}")

    layout = raw_layout()
    if layout:
        print(f"\nraw tree ({sum(layout.values())} frames):")
        for relative, count in layout.items():
            print(f"    {relative:<40s} {count:4d}")
    print("\nRun `cassa-doctor` for a full environment check.")


def raw_frames():
    """Every raw frame under ``RAW_DIR``, sorted, whatever the layout.

    The acquisition software delivers a night as a tree --
    ``<date>/BIAS|DARK|FLAT|LIGHT/<target>/`` -- while ``cassa-simulate`` writes
    one flat directory. Notebooks call this instead of globbing ``RAW_DIR/*.fits``
    so both read the same way, exactly as the pipeline itself reads them.
    """
    from cassa_photometry.paths import find_raw_frames

    return find_raw_frames(RAW_DIR)


def raw_layout():
    """``{directory relative to RAW_DIR: frame count}`` -- what is actually there."""
    from cassa_photometry.paths import raw_tree_summary

    if not os.path.isdir(RAW_DIR):
        return {}
    return raw_tree_summary(raw_frames(), RAW_DIR)


def require_dataset():
    """Fail early and clearly if the dataset has not been generated."""
    if raw_frames():
        return
    raise FileNotFoundError(
        f"No FITS frames in {RAW_DIR} or below it.\n\n"
        f"Either copy a night from the acquisition software in (its own\n"
        f"<date>/BIAS|DARK|FLAT|LIGHT/<target>/ tree is read as-is), or simulate\n"
        f"one -- a couple of minutes, and reproducible from a seed:\n\n"
        f"    cd {WORKSHOP_DIR}\n"
        f"    cassa-simulate --preset workshop --out .\n"
    )


def load_truth():
    """The truth catalogue and per-frame truth, as DataFrames."""
    import pandas as pd

    if not os.path.exists(TRUTH_SOURCES):
        raise FileNotFoundError(
            f"No truth file at {TRUTH_SOURCES}. It is written alongside the raw "
            f"frames by `cassa-simulate`."
        )
    return pd.read_csv(TRUTH_SOURCES), pd.read_csv(TRUTH_FRAMES)


def score_against_truth(catalog_csv, band, columns=("MAG_BEST", "MAG_APER", "MAG_ISO"),
                        bright_limit=16.5, match_arcsec=1.5):
    """Compare a catalog's magnitudes against the simulated truth.

    Returns a DataFrame with, per column, the median error, the scatter, and the
    **trend with brightness** -- which is the one that matters. A constant offset
    is absorbed by any zero point; a trend is not, and it is what distinguishes
    an isophotal magnitude from a total one.
    """
    import astropy.units as u
    import numpy as np
    import pandas as pd
    from astropy.coordinates import SkyCoord

    truth, _ = load_truth()
    stars = truth[truth["type"] == "star"]
    truth_coords = SkyCoord(stars["ra"].to_numpy() * u.deg, stars["dec"].to_numpy() * u.deg)

    catalog = pd.read_csv(catalog_csv)
    catalog_coords = SkyCoord(catalog["ALPHA_J2000"].to_numpy() * u.deg,
                              catalog["DELTA_J2000"].to_numpy() * u.deg)
    index, separation, _ = catalog_coords.match_to_catalog_sky(truth_coords)
    matched = np.asarray(separation < match_arcsec * u.arcsec)

    true_mag = stars[f"mag_{band}"].to_numpy()[index[matched]]
    rows = []
    for column in columns:
        if column not in catalog.columns:
            continue
        measured = catalog[column].to_numpy()[matched]
        good = np.isfinite(measured) & np.isfinite(true_mag) & (true_mag < bright_limit)
        if good.sum() < 5:
            continue
        error = measured[good] - true_mag[good]
        rows.append({
            "column": column,
            "n": int(good.sum()),
            "median_error": float(np.median(error)),
            "scatter": float(np.std(error)),
            "trend_mag_per_mag": float(np.polyfit(true_mag[good], error, 1)[0]),
        })
    return pd.DataFrame(rows)
