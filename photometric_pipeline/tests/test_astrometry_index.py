"""Selecting and fetching only the index files a field needs.

The point of this machinery is that a fresh clone can plate-solve without a
5 GB manual download, so the tests are about two properties: that the selection
is *correct* (never omits a file the solver needs), and that fetching is *safe*
(a truncated or corrupt download never masquerades as a usable index).
"""

import http.server
import os
import threading

import numpy as np
import pytest

from cassa_photometry.astrometry_index import (
    ALL_SKY,
    HttpIndexStore,
    IndexEntry,
    IndexManifest,
    LocalIndexStore,
    _cone_samples,
    _healpixes_in_cone,
    field_size_arcmin,
    healpix_available,
    required_indexes,
)
from cassa_photometry.config import load_config

NGC7331 = (339.266875, 34.415778)


def _entry(name, series, healpix, nside, lo, hi, size=1000):
    return IndexEntry(name=name, series=series, healpix=healpix, nside=nside,
                      scale_lo_arcmin=lo, scale_hi_arcmin=hi, bytes=size,
                      url_path=name)


@pytest.fixture
def manifest():
    """A small stand-in with the same structure as the shipped one."""
    return IndexManifest([
        _entry("index-4202-14.fits", 4202, 14, 2, 4.0, 5.6, 109_000_000),
        _entry("index-4203-14.fits", 4203, 14, 2, 5.6, 8.0, 55_000_000),
        _entry("index-4203-15.fits", 4203, 15, 2, 5.6, 8.0, 55_000_000),
        _entry("index-4204-14.fits", 4204, 14, 2, 8.0, 11.0, 27_000_000),
        _entry("index-4205-03.fits", 4205, 3, 1, 11.0, 16.0, 55_000_000),
        _entry("index-4212.fits", 4212, ALL_SKY, 0, 120.0, 170.0, 1_000_000),
    ])


# --- Field geometry -----------------------------------------------------------

def test_field_size_uses_the_short_axis():
    """Conservative on purpose: never ask for quads the field cannot contain."""
    assert field_size_arcmin(0.598, 2048, 1400) == pytest.approx(1400 * 0.598 / 60)
    assert field_size_arcmin(0.598, 1400, 2048) == pytest.approx(1400 * 0.598 / 60)


def test_field_size_admits_when_it_cannot_say():
    assert field_size_arcmin(None, 1024, 1024) is None
    assert field_size_arcmin(0.0, 1024, 1024) is None
    assert field_size_arcmin(0.598, 0, 1024) is None


# --- Selection ----------------------------------------------------------------

def test_only_indexes_whose_quads_fit_the_field_are_selected(manifest):
    """A 10' field can use 5.6-11' quads, not 120-170' ones."""
    selected = required_indexes(*NGC7331, 0.598, 1024, 1024, manifest=manifest,
                                scale_lo_frac=0.30, scale_hi_frac=1.0)
    names = {e.name for e in selected}
    assert "index-4203-14.fits" in names
    assert "index-4204-14.fits" in names
    assert "index-4212.fits" not in names, "an all-sky index of the wrong scale"


def test_a_bigger_field_reaches_coarser_indexes(manifest):
    small = required_indexes(*NGC7331, 0.598, 768, 768, manifest=manifest,
                             scale_lo_frac=0.30)
    large = required_indexes(*NGC7331, 0.598, 4096, 4096, manifest=manifest,
                             scale_lo_frac=0.30)
    assert max(e.scale_hi_arcmin for e in large) > max(e.scale_hi_arcmin for e in small)


def test_only_the_tiles_covering_the_pointing_are_selected(manifest):
    """The whole saving: a field needs its own sky tile, not all 48."""
    if not healpix_available():
        pytest.skip("HEALPix bindings unavailable; selection falls back to a superset")
    selected = required_indexes(*NGC7331, 0.598, 1024, 1024, manifest=manifest,
                                scale_lo_frac=0.30)
    tiles = {e.healpix for e in selected if e.healpix >= 0}
    assert tiles == {14}


def test_an_all_sky_index_is_kept_whenever_its_scale_matches(manifest):
    selected = required_indexes(*NGC7331, 60.0, 1024, 1024, manifest=manifest,
                                scale_lo_frac=0.10, scale_hi_frac=1.0)
    assert "index-4212.fits" in {e.name for e in selected}


def test_selection_is_empty_without_a_pointing_or_a_scale(manifest):
    """Not an error: a frame with no pointing is a normal thing to meet, and the
    caller falls back to whatever local set exists."""
    assert required_indexes(None, None, 0.598, 1024, 1024, manifest=manifest) == []
    assert required_indexes(*NGC7331, None, 1024, 1024, manifest=manifest) == []


def test_finest_quads_come_first(manifest):
    selected = required_indexes(*NGC7331, 0.598, 2048, 2048, manifest=manifest,
                                scale_lo_frac=0.10)
    scales = [e.scale_lo_arcmin for e in selected]
    assert scales == sorted(scales)


# --- HEALPix ------------------------------------------------------------------

def test_a_cone_is_sampled_on_its_boundary_not_just_its_centre():
    samples = _cone_samples(10.0, 20.0, 3.0, n_boundary=36)
    assert len(samples) == 37
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    centre = SkyCoord(10.0 * u.deg, 20.0 * u.deg)
    edge = SkyCoord([s[0] for s in samples[1:]] * u.deg,
                    [s[1] for s in samples[1:]] * u.deg)
    assert np.allclose(centre.separation(edge).deg, 3.0, atol=1e-6)


def test_a_pointing_near_a_tile_boundary_pulls_in_its_neighbour():
    """Centre-only selection silently misses these, and the solve then fails."""
    if not healpix_available():
        pytest.skip("HEALPix bindings unavailable")
    rng = np.random.default_rng(0)
    straddling = 0
    for _ in range(200):
        ra = rng.uniform(0, 360)
        dec = np.degrees(np.arcsin(rng.uniform(-1, 1)))
        if len(_healpixes_in_cone(ra, dec, 3.0, 2)) > 1:
            straddling += 1
    assert straddling > 10, "boundary sampling never found a multi-tile pointing"


def test_without_healpix_bindings_selection_is_a_superset(manifest, monkeypatch):
    """The fallback may download more, never less: the solve stays correct."""
    import cassa_photometry.astrometry_index as module

    exact = required_indexes(*NGC7331, 0.598, 1024, 1024, manifest=manifest,
                             scale_lo_frac=0.10)
    monkeypatch.setattr(module, "_healpix_of", lambda *a, **k: None)
    fallback = required_indexes(*NGC7331, 0.598, 1024, 1024, manifest=manifest,
                                scale_lo_frac=0.10)
    assert {e.name for e in exact} <= {e.name for e in fallback}


# --- The shipped manifest -----------------------------------------------------

def test_the_shipped_manifest_describes_the_public_index_server():
    manifest = IndexManifest.default()
    assert len(manifest) > 250, "manifest missing or truncated"
    series = {e.series for e in manifest.entries}
    assert {4203, 4204, 4205, 4206} <= series
    # nside follows from the file count per series; getting it wrong would send
    # the selector to the wrong tile.
    by_name = {e.name: e for e in manifest.entries}
    assert by_name["index-4203-14.fits"].nside == 2
    assert by_name["index-4205-03.fits"].nside == 1
    assert all(e.bytes > 0 for e in manifest.entries)


def test_a_cassa_field_needs_a_small_fraction_of_the_set():
    """The claim the whole feature rests on.

    The exact figure depends on whether the native tile numbering is available:
    a conda install fetches ~246 MB, a pip-only one ~543 MB via the geometric
    fallback. Either is a small fraction of the ~34 GB on the server, which is
    what matters.
    """
    manifest = IndexManifest.default()
    selected = required_indexes(*NGC7331, 0.598, 2048, 1400, manifest=manifest,
                                scale_lo_frac=0.30)
    needed = sum(e.bytes for e in selected)
    total = sum(e.bytes for e in manifest.entries)
    assert 0 < needed / total < 0.03
    assert needed / 1e6 < 700, "first-pass download should stay well under a GB"


def test_the_geometric_fallback_is_tight_as_well_as_correct(monkeypatch):
    """A pip-only install has no astrometry.net bindings, and no other HEALPix
    library reproduces their tile numbering. The manifest therefore stores each
    tile's sky position, so the fallback selects neighbours rather than
    everything."""
    import cassa_photometry.astrometry_index as module

    manifest = IndexManifest.default()
    exact = required_indexes(*NGC7331, 0.598, 2048, 1400, manifest=manifest,
                             scale_lo_frac=0.30)

    monkeypatch.setattr(module, "_healpix_of", lambda *a, **k: None)
    fallback = required_indexes(*NGC7331, 0.598, 2048, 1400, manifest=manifest,
                                scale_lo_frac=0.30)

    assert {e.name for e in exact} <= {e.name for e in fallback}, "not a superset"
    # Tight: a handful of neighbouring tiles, not all 48 of every series.
    assert len(fallback) < 4 * len(exact)


def test_the_manifest_records_where_each_tile_is():
    manifest = IndexManifest.default()
    tiled = [e for e in manifest.entries if not e.is_all_sky]
    assert all(e.center_ra is not None and e.radius_deg for e in tiled)

    tile = next(e for e in tiled if e.name == "index-4203-14.fits")
    assert tile.reaches(*NGC7331, 3.0)
    assert not tile.reaches(0.0, -80.0, 3.0)


def test_manifest_round_trips(tmp_path, manifest):
    path = tmp_path / "m.json"
    manifest.save(str(path))
    again = IndexManifest.load(str(path))
    assert [e.name for e in again.entries] == [e.name for e in manifest.entries]
    assert again.entries[0].scale_lo_arcmin == manifest.entries[0].scale_lo_arcmin


# --- Stores -------------------------------------------------------------------

def test_local_store_returns_only_what_exists(tmp_path, manifest):
    (tmp_path / "index-4203-14.fits").write_bytes(b"x")
    store = LocalIndexStore(str(tmp_path))
    paths = store.ensure(manifest.entries)
    assert len(paths) == 1 and paths[0].endswith("index-4203-14.fits")


def test_an_existing_local_set_is_preferred_over_downloading(tmp_path, manifest):
    local = tmp_path / "local"
    local.mkdir()
    (local / "index-4203-14.fits").write_bytes(b"x")
    store = HttpIndexStore("http://127.0.0.1:1/", str(tmp_path / "cache"),
                           local_dirs=[str(local)])
    assert store.path_of(manifest.entries[1]).startswith(str(local))


class _QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):  # keep the test output clean
        pass


@pytest.fixture
def index_server(tmp_path):
    """A real HTTP server over a directory of fake index files."""
    root = tmp_path / "server"
    root.mkdir()
    (root / "index-4203-14.fits").write_bytes(b"A" * 4096)

    handler = type("H", (_QuietHandler,), {"directory": str(root)})
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0),
                                             lambda *a, **k: handler(*a, directory=str(root), **k))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_port}/", root
    server.shutdown()


def test_a_missing_file_is_downloaded_and_cached(tmp_path, index_server, manifest):
    base_url, _ = index_server
    cache = tmp_path / "cache"
    store = HttpIndexStore(base_url, str(cache))
    entry = manifest.entries[1]

    paths = store.ensure([entry])
    assert len(paths) == 1
    assert os.path.getsize(paths[0]) == 4096
    # ...and a second call uses the cache rather than the network.
    store.base_url = "http://127.0.0.1:1/"
    assert store.ensure([entry]) == paths


def test_a_checksum_mismatch_is_rejected_rather_than_used(tmp_path, index_server, manifest):
    """A corrupt index makes a solve fail in a way nothing else would explain,
    so it must never be left in the cache looking valid."""
    base_url, _ = index_server
    cache = tmp_path / "cache"
    entry = manifest.entries[1]
    entry.sha256 = "0" * 64
    store = HttpIndexStore(base_url, str(cache))

    assert store.ensure([entry]) == []
    assert not os.path.exists(os.path.join(str(cache), entry.name))
    assert not os.path.exists(os.path.join(str(cache), entry.name + ".part"))


def test_a_missing_download_does_not_stop_the_run(tmp_path, manifest):
    store = HttpIndexStore("http://127.0.0.1:1/", str(tmp_path / "cache"))
    assert store.ensure([manifest.entries[0]]) == []


def test_downloads_can_be_disabled(tmp_path, manifest, pipeline_logs):
    store = HttpIndexStore("http://127.0.0.1:1/", str(tmp_path / "cache"), download=False)
    assert store.ensure([manifest.entries[0]]) == []
    assert "cassa-index-fetch" in pipeline_logs.text


def test_eviction_never_removes_a_file_this_run_needs(tmp_path, manifest):
    cache = tmp_path / "cache"
    cache.mkdir()
    keep = manifest.entries[1]
    (cache / keep.name).write_bytes(b"A" * 2048)
    (cache / "index-9999-00.fits").write_bytes(b"B" * 2048)

    store = HttpIndexStore("http://127.0.0.1:1/", str(cache), cache_gb=1e-9)
    store._evict_for(0, keep=[keep])

    assert (cache / keep.name).exists(), "evicted a file the run needs"
    assert not (cache / "index-9999-00.fits").exists()


def test_build_store_prefers_the_configured_local_directory(tmp_path, monkeypatch):
    monkeypatch.setenv("CASSA_ASTROMETRY_INDEX", str(tmp_path))
    from cassa_photometry.astrometry_index import build_store

    store = build_store(load_config())
    assert str(tmp_path) in store.local_dirs
