"""The on-disk layout shared by every phase.

A run lives in one *work directory* holding a directory per phase::

    work/
      phase1/   calibrated frames        (cassa-calibrate)
      phase2/   master stacks + WCS      (cassa-integrate)
      phase3/   flux-calibrated + catalogs (cassa-photometry)
      phase4/   diagnostics report       (cassa-diagnose)

Each phase writes to its own directory and reads from the previous one, so a
re-run replaces that phase's products in place instead of accumulating copies.
"""

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
