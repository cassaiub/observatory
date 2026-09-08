"""The on-disk layout shared by every phase.

A run lives in one *work directory* holding a directory per phase::

    work/
      phase1/   calibrated frames        (cassa-calibrate)
      phase2/   master stacks + WCS      (cassa-integrate)
      phase3/   flux-calibrated + catalogs (cassa-photometry)
      phase4/   diagnostics report       (cassa-diagnose)

Each phase writes to its own directory and reads from the previous one, so a
re-run replaces that phase's products in place instead of accumulating copies.

The *raw* tree is the one layout the pipeline does not own -- the acquisition
software chooses it -- so :func:`find_raw_frames` reads it whatever shape it
arrives in, and every phase that takes a raw directory goes through it.
"""

import glob
import os

PHASE_NAMES = {1: "phase1", 2: "phase2", 3: "phase3", 4: "phase4"}


def phase_dir(work_dir, phase):
    """Return ``<work_dir>/phaseN``."""
    return os.path.join(os.path.abspath(work_dir), PHASE_NAMES[phase])


def sibling_phase_dir(other_phase_dir, phase):
    """Return the ``phaseN`` directory next to an existing phase directory.

    ``sibling_phase_dir("work/phase1", 2)`` -> ``work/phase2``.
    """
    return phase_dir(os.path.dirname(os.path.abspath(other_phase_dir)), phase)


def find_phase_dir(other_phase_dir, phase):
    """Sibling ``phaseN`` directory if it exists, else ``None``."""
    candidate = sibling_phase_dir(other_phase_dir, phase)
    return candidate if os.path.isdir(candidate) else None


#: Filename glob matching the FITS spellings acquisition software writes
#: (``.fits``, ``.fit``, ``.fts``, and their ``.gz`` forms).
RAW_FITS_GLOB = "*.f*t*"


def find_raw_frames(root, pattern=RAW_FITS_GLOB):
    """Every raw frame at or below ``root``, sorted, whatever the layout.

    The acquisition software delivers a night as a tree::

        raw/
          20260903/
            BIAS/untargeted/  20260903T054529.715412_..._BIAS_....fits
            DARK/untargeted/  ...
            FLAT/untargeted/  ...
            LIGHT/m22/        20260903T043707.171323_m22_..._LIGHT_....fits

    while ``cassa-simulate`` and hand-assembled sets put every frame in one flat
    directory. Both are read by walking the tree, so no caller has to know which
    one it was handed, and a raw argument that names a single file is returned as
    itself.

    The directory names are *not* what classifies a frame -- the instrument
    profile reads ``IMAGETYP`` for that, as it always has. The tree only says
    where the files are, and a night sorted into the wrong folder is still
    reduced as whatever its header says it is.
    """
    if os.path.isfile(root):
        return [os.path.abspath(root)]
    return sorted(
        p for p in glob.glob(os.path.join(root, "**", pattern), recursive=True)
        if os.path.isfile(p)
    )


def raw_tree_summary(paths, root):
    """``{directory relative to root: frame count}``, in path order.

    For logging what a raw argument actually resolved to: a flat directory
    prints one line, a night from the acquisition software prints one per
    ``<date>/<type>/<target>`` leaf.
    """
    counts = {}
    root = os.path.abspath(root)
    for path in paths:
        directory = os.path.dirname(os.path.abspath(path))
        try:
            key = os.path.relpath(directory, root)
        except ValueError:            # different drive on Windows
            key = directory
        counts[key] = counts.get(key, 0) + 1
    return {key: counts[key] for key in sorted(counts)}
