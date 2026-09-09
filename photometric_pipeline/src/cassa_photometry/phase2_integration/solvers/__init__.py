"""Plate-solving backends.

``get_solver`` resolves ``phase2.solver`` to something that can actually run
here, so a machine with no ``solve-field`` binary still solves and a machine
with one still uses it.
"""

from cassa_photometry.phase2_integration.solvers.astap import AstapSolver
from cassa_photometry.phase2_integration.solvers.base import (
    SolveHints,
    Solver,
    SolveResult,
    solved_pixel_scale,
)
from cassa_photometry.phase2_integration.solvers.inprocess import InProcessSolver
from cassa_photometry.phase2_integration.solvers.solvefield import SolveFieldSolver

#: Config name -> backend, in the order ``auto`` prefers them.
#:
#: ASTAP leads because it installs everywhere this pipeline runs (including ARM
#: Linux, where the other two have no build at all), needs no Python package,
#: and fetches only the sky tiles a field uses. ``solve-field`` follows because
#: it is the reference implementation and the only backend that supplies its own
#: matched-star table. The in-process solver is last: it is what keeps a
#: pip-only machine working when neither binary can be installed.
BACKENDS = {
    "astap": AstapSolver,
    "solve-field": SolveFieldSolver,
    "astrometry-py": InProcessSolver,
}

__all__ = [
    "SolveHints", "SolveResult", "Solver", "solved_pixel_scale",
    "SolveFieldSolver", "InProcessSolver", "AstapSolver", "BACKENDS",
    "get_solver", "available_backends",
]


def _can_run(cls, config):
    """Whether a backend can run, passing config to those that need it.

    ASTAP's availability depends on ``phase2.astap_path`` as well as PATH, so it
    accepts a config; the other two do not. Handling the difference here keeps
    the Solver interface uniform rather than forcing every backend to take an
    argument it ignores.
    """
    try:
        return cls.available(config)
    except TypeError:
        return cls.available()


def available_backends(config=None):
    """Names of the backends that can run on this machine."""
    return [name for name, cls in BACKENDS.items() if _can_run(cls, config)]


def get_solver(logger, config):
    """The solver this configuration asks for, falling back when it cannot run."""
    requested = (getattr(config.phase2, "solver", "auto") or "auto").strip().lower()

    if requested == "auto":
        for cls in BACKENDS.values():
            if _can_run(cls, config):
                return cls(logger, config)
        logger.error(
            "No plate solver is available. Install ASTAP "
            "(apt install astap-cli, or https://www.hnsky.org/astap.htm), "
            "Astrometry.net (conda install -c conda-forge astrometry), or the "
            'in-process solver (pip install "cassa-photometry[solver]"). '
            "Run cassa-doctor for details."
        )
        return None

    if requested not in BACKENDS:
        raise KeyError(
            f"Unknown solver {requested!r}. Available: auto, {', '.join(BACKENDS)}."
        )

    cls = BACKENDS[requested]
    if not _can_run(cls, config):
        logger.warning(
            "Solver %r was requested but is not installed on this machine.", requested
        )
        for name, other in BACKENDS.items():
            if _can_run(other, config):
                logger.warning("Falling back to %r.", name)
                return other(logger, config)
        return None
    return cls(logger, config)
