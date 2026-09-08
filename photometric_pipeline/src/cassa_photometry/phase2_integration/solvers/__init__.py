"""Plate-solving backends.

``get_solver`` resolves ``phase2.solver`` to something that can actually run
here, so a machine with no ``solve-field`` binary still solves and a machine
with one still uses it.
"""

from cassa_photometry.phase2_integration.solvers.base import (
    SolveHints,
    Solver,
    SolveResult,
    solved_pixel_scale,
)
from cassa_photometry.phase2_integration.solvers.inprocess import InProcessSolver
from cassa_photometry.phase2_integration.solvers.solvefield import SolveFieldSolver

#: Config name -> backend, in the order ``auto`` prefers them.
BACKENDS = {
    "solve-field": SolveFieldSolver,
    "astrometry-py": InProcessSolver,
}

__all__ = [
    "SolveHints", "SolveResult", "Solver", "solved_pixel_scale",
    "SolveFieldSolver", "InProcessSolver", "BACKENDS",
    "get_solver", "available_backends",
]


def available_backends():
    """Names of the backends that can run on this machine."""
    return [name for name, cls in BACKENDS.items() if cls.available()]


def get_solver(logger, config):
    """The solver this configuration asks for, falling back when it cannot run."""
    requested = (getattr(config.phase2, "solver", "auto") or "auto").strip().lower()

    if requested == "auto":
        for cls in BACKENDS.values():
            if cls.available():
                return cls(logger, config)
        logger.error(
            "No plate solver is available. Install Astrometry.net "
            '(conda install -c conda-forge astrometry) or the in-process solver '
            '(pip install "cassa-photometry[solver]"). Run cassa-doctor for details.'
        )
        return None

    if requested not in BACKENDS:
        raise KeyError(
            f"Unknown solver {requested!r}. Available: auto, {', '.join(BACKENDS)}."
        )

    cls = BACKENDS[requested]
    if not cls.available():
        logger.warning(
            "Solver %r was requested but is not installed on this machine.", requested
        )
        for name, other in BACKENDS.items():
            if other.available():
                logger.warning("Falling back to %r.", name)
                return other(logger, config)
        return None
    return cls(logger, config)
