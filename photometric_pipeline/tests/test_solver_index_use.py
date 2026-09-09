"""Only backends that read index files should cause index files to be fetched.

Phase 2 used to select and download Astrometry.net index files on every solve,
whatever the backend. ASTAP ignores them -- it carries its own star database --
so an ASTAP user paid ~246 MB per field for data handed to a solver that never
opens it, on top of the ~6 MB of tiles it does use. Reported from a real run::

    -> Solving with 8 index file(s); SolveHints(...)
    -> ASTAP database: 10 tile(s) already cached.

Both lines for one solve, and only the second one mattered.
"""

from unittest import mock

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.phase2_integration.solvers.astap import AstapSolver
from cassa_photometry.phase2_integration.solvers.base import Solver
from cassa_photometry.phase2_integration.solvers.inprocess import InProcessSolver
from cassa_photometry.phase2_integration.solvers.solvefield import SolveFieldSolver
from cassa_photometry.phase2_integration.wcs import WCSSolver


def test_backends_declare_whether_they_read_index_files():
    assert Solver.uses_index_files is True, "the safe default is to fetch"
    assert AstapSolver.uses_index_files is False
    # These two genuinely solve against index files.
    assert SolveFieldSolver.uses_index_files is True
    assert InProcessSolver.uses_index_files is True


def _wcs_solver_with(backend):
    with mock.patch("cassa_photometry.phase2_integration.wcs.get_solver",
                    return_value=backend):
        return WCSSolver(get_logger("t"), load_config())


def test_astap_selects_no_index_files_at_all():
    """The saving is the point: nothing should even be *selected*, so no
    manifest lookup and no download can follow."""
    solver = _wcs_solver_with(AstapSolver(get_logger("t"), load_config()))

    with mock.patch("cassa_photometry.phase2_integration.wcs.required_indexes") as chooser:
        assert solver._indexes_for(mock.Mock(), 0.3) == []
    chooser.assert_not_called()


def test_a_backend_that_uses_indexes_still_gets_them():
    solver = _wcs_solver_with(SolveFieldSolver(get_logger("t"), load_config()))

    with mock.patch("cassa_photometry.phase2_integration.wcs.required_indexes",
                    return_value=[]) as chooser:
        solver._indexes_for(mock.Mock(), 0.3)
    chooser.assert_called_once()


def test_astap_is_not_warned_about_a_missing_index_directory(tmp_path, caplog):
    """"WCS solving will fail" is wrong and expensive advice for an ASTAP user:
    it sends them to download several gigabytes they do not need."""
    config = load_config()
    config.phase2.astrometry_index_dir = str(tmp_path / "nothing-here")

    with caplog.at_level("WARNING"):
        _wcs_solver_with(AstapSolver(get_logger("t"), config))

    assert "cassa-index-fetch" not in caplog.text
