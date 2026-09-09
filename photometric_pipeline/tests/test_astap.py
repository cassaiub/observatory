"""The ASTAP backend, and the on-demand tile store behind it.

Two properties matter more than any individual assertion here:

1. **Both install paths produce the same products.** ASTAP supplies no matched
   star table, so ``ASTRMS`` -- which phase 3 sizes its cross-match radius from
   -- has to be measured by the pipeline instead. If that stopped working, a run
   would silently fall back to a fixed 2" radius depending on which solver
   happened to be installed.
2. **A field's tiles are worked out, not downloaded wholesale.** The published
   archive is 859 MB; a field needs about 6 MB of it. The selection arithmetic
   is what makes that difference, so it is tested directly.

Nothing here needs ASTAP installed or a network: the binary is faked where the
subprocess boundary is what is under test, and the tile maths is pure.
"""

import os

import numpy as np
import pytest

from cassa_photometry import astap_db
from cassa_photometry.config import load_config
from cassa_photometry.phase2_integration.solvers import BACKENDS, get_solver
from cassa_photometry.phase2_integration.solvers.astap import AstapSolver, find_binary

# --- Tile selection -----------------------------------------------------------

def test_the_sky_is_covered_exactly_once():
    """36 bands of 5 degrees, and the per-band tile counts are the real ones.

    If this table drifts, every selection is wrong by a tile and fields near the
    edges silently lose coverage.
    """
    assert len(astap_db.TILES_PER_BAND) == 36
    assert astap_db.BAND_HEIGHT_DEG == 5.0
    # Symmetric about the equator: the sky does not care which pole you start at.
    assert list(astap_db.TILES_PER_BAND) == list(reversed(astap_db.TILES_PER_BAND))
    # Poles hold one tile; the equator holds the most.
    assert astap_db.TILES_PER_BAND[0] == 1
    assert max(astap_db.TILES_PER_BAND) == 69


@pytest.mark.parametrize(("dec", "band"), [
    (-90.0, 1), (-89.9, 1), (0.0, 19), (34.4158, 25), (89.9, 36), (90.0, 36),
])
def test_declination_maps_to_the_right_band(dec, band):
    assert astap_db.band_of(dec) == band


def test_right_ascension_wraps():
    """359.9 and 0.1 degrees are neighbours, not opposite ends of the sky."""
    band = 19
    assert astap_db.tile_of(0.0, band) == 1
    assert astap_db.tile_of(360.0, band) == 1
    assert astap_db.tile_of(359.999, band) == astap_db.TILES_PER_BAND[band - 1]


def test_the_selection_matches_the_field_that_was_measured():
    """The exact tiles a real NGC7331 solve needed -- 10 of 1476."""
    tiles = astap_db.required_tiles(339.2669, 34.4158, "d50")
    assert tiles == [
        "d50_0101", "d50_2456", "d50_2457", "d50_2458",
        "d50_2552", "d50_2553", "d50_2554",
        "d50_2649", "d50_2650", "d50_2651",
    ]


def test_the_marker_tile_is_always_included():
    """Without <series>_0101 ASTAP reports "no star database found" however many
    other tiles are present -- it is how the program recognises a database."""
    for dec in (-80.0, 0.0, 12.5, 80.0):
        assert "d50_0101" in astap_db.required_tiles(10.0, dec, "d50")


def test_a_pole_field_does_not_ask_for_bands_off_the_sky():
    tiles = astap_db.required_tiles(0.0, 89.0, "d50")
    bands = {int(t.split("_")[1][:2]) for t in tiles}
    assert max(bands) <= 36 and min(bands) >= 1


def test_no_pointing_selects_nothing_rather_than_guessing():
    assert astap_db.required_tiles(None, None) == []


@pytest.mark.parametrize(("fov", "series"), [
    (0.13, "d50"),   # a CASSA 14-inch field
    (0.36, "d50"),   # the CASSA 8-inch
    (12.0, "g05"),   # a wide-field rig
])
def test_series_is_chosen_by_field_size(fov, series):
    assert astap_db.series_for_fov(fov) == series


# --- The local store ----------------------------------------------------------

def test_tiles_already_on_disk_are_not_refetched(tmp_path, monkeypatch):
    """Repointing must cost only the tiles the new field adds."""
    for stem in astap_db.required_tiles(339.2669, 34.4158, "d50"):
        (tmp_path / f"{stem}.1476").write_bytes(b"x")

    store = astap_db.AstapTileStore(str(tmp_path))

    def explode(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("a cached field must not trigger a download")

    monkeypatch.setattr(store, "_fetch", explode)
    assert store.ensure(339.2669, 34.4158, "d50") == str(tmp_path)


def test_missing_tiles_are_reported_when_downloading_is_disabled(tmp_path):
    """The air-gapped case: say which tiles are needed rather than hanging."""
    store = astap_db.AstapTileStore(str(tmp_path), download=False)
    with pytest.raises(astap_db.AstapDatabaseError, match="missing"):
        store.ensure(339.2669, 34.4158, "d50")


def test_an_unfetchable_series_is_named(tmp_path):
    """d80 is published only as a .deb and a Windows installer, neither of which
    supports per-member range reads. Saying so beats a confusing HTTP error."""
    store = astap_db.AstapTileStore(str(tmp_path))
    with pytest.raises(astap_db.AstapDatabaseError, match="d80"):
        store._archive_url("d80")


def test_the_cache_directory_resolves_like_the_index_cache(monkeypatch):
    config = load_config()
    monkeypatch.setenv("CASSA_ASTAP_DB", "/tmp/astap-tiles")
    assert astap_db.resolve_db_dir(config) == "/tmp/astap-tiles"

    config.phase2.astap_db_dir = "/explicit/path"
    assert astap_db.resolve_db_dir(config) == "/explicit/path"


# --- The backend --------------------------------------------------------------

def test_astap_is_registered_and_preferred():
    """`auto` should reach for ASTAP first: it is the only backend that installs
    on every platform this pipeline supports."""
    assert "astap" in BACKENDS
    assert list(BACKENDS) == ["astap", "solve-field", "astrometry-py"]


def test_an_absent_binary_is_not_available():
    assert find_binary("/nonexistent/astap_cli") is None
    assert not AstapSolver.available(_config_with_path("/nonexistent/astap_cli"))


def _config_with_path(path):
    config = load_config()
    config.phase2.astap_path = path
    return config


def test_hints_are_converted_to_astap_conventions():
    """RA in HOURS and declination as south-polar distance. Getting either wrong
    produces a failed solve with no hint as to why."""
    from cassa_photometry.logging_utils import get_logger
    from cassa_photometry.phase2_integration.solvers.base import SolveHints

    solver = AstapSolver(get_logger("t"), load_config())
    args = solver._hint_args(SolveHints(
        ra_deg=339.2669, dec_deg=34.4158, pixel_scale=0.591,
        naxis1=1024, naxis2=1024, radius_deg=10.0))

    assert args[args.index("-ra") + 1] == f"{339.2669 / 15.0:.6f}"
    assert args[args.index("-spd") + 1] == f"{34.4158 + 90.0:.6f}"
    # Field HEIGHT in degrees, from the frame's own scale.
    assert float(args[args.index("-fov") + 1]) == pytest.approx(
        1024 * 0.591 / 3600.0, rel=1e-4)


def test_missing_coverage_is_distinguished_from_an_unsolvable_frame():
    """ASTAP says "Error reading star database" for the first and "no solution"
    for the second. Collapsing them would send a user hunting for a database
    they already have."""
    from cassa_photometry.logging_utils import get_logger

    solver = AstapSolver(get_logger("t"), load_config())
    completed = type("P", (), {"stdout": "No solution found!", "stderr": ""})()

    coverage = solver._failure_message(
        {"ERROR": "Error reading star database."}, completed)
    assert "no star database covering this field" in coverage

    unsolvable = solver._failure_message({}, completed)
    assert "No solution found" in unsolvable


def test_a_solved_ini_becomes_a_usable_wcs(tmp_path):
    """The .ini is ASTAP's real output; the .wcs file embeds the whole original
    header (CONTINUE cards included), which astropy refuses to re-serialise."""
    from astropy.io import fits
    from astropy.wcs import WCS

    from cassa_photometry.logging_utils import get_logger

    frame = tmp_path / "m.fits"
    fits.PrimaryHDU(np.zeros((16, 16), dtype=np.float32)).writeto(frame)

    solver = AstapSolver(get_logger("t"), load_config())
    header = solver._header({
        "CRPIX1": "5.125E+002", "CRPIX2": "5.125E+002",
        "CRVAL1": "3.3917396788E+002", "CRVAL2": "3.4411151984E+001",
        "CDELT1": "-1.6404563398E-004", "CDELT2": "1.6395231995E-004",
        "CD1_1": "-1.6404561805E-004", "CD1_2": "7.2287359895E-008",
        "CD2_1": "-4.7840110074E-009", "CD2_2": "1.6395231988E-004",
    }, str(frame))

    assert header is not None
    wcs = WCS(header)
    assert wcs.has_celestial
    from astropy.wcs.utils import proj_plane_pixel_scales
    scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial)) * 3600)
    assert scale == pytest.approx(0.5905, abs=0.002)


def test_an_ini_without_a_solution_yields_no_header(tmp_path):
    from astropy.io import fits

    from cassa_photometry.logging_utils import get_logger

    frame = tmp_path / "m.fits"
    fits.PrimaryHDU(np.zeros((4, 4), dtype=np.float32)).writeto(frame)
    solver = AstapSolver(get_logger("t"), load_config())
    assert solver._header({"PLTSOLVD": "F"}, str(frame)) is None


def test_get_solver_falls_through_when_astap_is_absent():
    """The install chain is ASTAP -> solve-field -> in-process; asking for a
    backend that is not installed must fall back rather than fail."""
    from cassa_photometry.logging_utils import get_logger

    config = load_config()
    config.phase2.solver = "astap"
    config.phase2.astap_path = "/nonexistent/astap_cli"
    # Whatever is installed here, this must not raise.
    get_solver(get_logger("t"), config)


def test_the_configured_binary_path_is_honoured(tmp_path):
    fake = tmp_path / "astap_cli"
    fake.write_text("#!/bin/sh\nexit 0\n")
    os.chmod(fake, 0o755)
    assert find_binary(str(fake)) == str(fake)
    assert AstapSolver.available(_config_with_path(str(fake)))
