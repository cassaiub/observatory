"""Fetch only the astrometry index files a field actually needs.

A full Astrometry.net index set is about 5 GB, and downloading it by hand is the
single largest barrier to anyone running this pipeline on their own machine. It
is also almost entirely waste: for a given field the solver reads a handful of
those files and ignores the rest. Measured on the workshop field, a CASSA 8-inch
pointing needs roughly 3% of the set.

So the pipeline works out which files a frame requires and fetches those. The
files are static public data served over plain HTTPS from
``data.astrometry.net``, so nothing has to be hosted and no account is needed --
a fresh clone runs with no configuration.

Selection is exact rather than heuristic:

* **Which scales.** Astrometry.net matches four-star "quads" whose size must sit
  between roughly 10% and 100% of the field. Each index file states the quad
  scale range it holds, so the ones that overlap the window are the ones needed.
* **Which sky tiles.** The tiled series divide the sky into HEALPix cells. The
  cells a search cone touches are found by sampling the cone's boundary, which
  correctly picks up the extra tiles a pointing near a boundary needs.

The pipeline's own ``resolve_astrometry_index_dir()`` still wins when it holds a
populated set, so an existing installation is unaffected.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import urllib.error
import urllib.request

import numpy as np

from cassa_photometry.logging_utils import get_logger

#: Shipped manifest describing every index file the default server offers.
MANIFEST_PATH = os.path.join(os.path.dirname(__file__), "data", "index_manifest.json")

#: HEALPix id used by an all-sky index file, which covers every pointing.
ALL_SKY = -1

#: Download chunk size, and how many times a failed transfer is retried.
CHUNK_BYTES = 1 << 20
MAX_RETRIES = 4


class IndexEntry:
    """One index file: where it sits on the sky, and what quad scales it holds."""

    __slots__ = ("name", "series", "healpix", "nside", "scale_lo_arcmin",
                 "scale_hi_arcmin", "bytes", "sha256", "url_path",
                 "center_ra", "center_dec", "radius_deg")

    def __init__(self, name, series, healpix, nside, scale_lo_arcmin,
                 scale_hi_arcmin, bytes=0, sha256=None, url_path=None,
                 center_ra=None, center_dec=None, radius_deg=None):
        self.name = name
        self.series = int(series)
        self.healpix = int(healpix)
        self.nside = int(nside)
        self.scale_lo_arcmin = float(scale_lo_arcmin)
        self.scale_hi_arcmin = float(scale_hi_arcmin)
        self.bytes = int(bytes or 0)
        self.sha256 = sha256
        #: Path under the base URL, e.g. "4200/index-4203-14.fits".
        self.url_path = url_path or name
        #: Where this tile sits on the sky, so selection works without the
        #: astrometry.net bindings (which use their own tile numbering that
        #: nothing else reproduces).
        self.center_ra = center_ra
        self.center_dec = center_dec
        self.radius_deg = radius_deg

    @classmethod
    def from_dict(cls, data):
        return cls(
            name=data["name"], series=data["series"], healpix=data["healpix"],
            nside=data.get("nside", 1), scale_lo_arcmin=data["scale_lo_arcmin"],
            scale_hi_arcmin=data["scale_hi_arcmin"], bytes=data.get("bytes", 0),
            sha256=data.get("sha256"), url_path=data.get("url_path"),
            center_ra=data.get("center_ra"), center_dec=data.get("center_dec"),
            radius_deg=data.get("radius_deg"),
        )

    def to_dict(self):
        return {
            "name": self.name, "series": self.series, "healpix": self.healpix,
            "nside": self.nside, "scale_lo_arcmin": self.scale_lo_arcmin,
            "scale_hi_arcmin": self.scale_hi_arcmin, "bytes": self.bytes,
            "sha256": self.sha256, "url_path": self.url_path,
            "center_ra": self.center_ra, "center_dec": self.center_dec,
            "radius_deg": self.radius_deg,
        }

    @property
    def is_all_sky(self):
        return self.healpix == ALL_SKY

    def reaches(self, ra_deg, dec_deg, radius_deg):
        """Whether a search cone can touch this tile.

        Used when the native tile numbering is unavailable. Comparing centre
        separation against the sum of the radii is a *superset* test -- a tile
        may be included that the cone does not truly overlap -- so the solve is
        never wrong, only the download slightly larger. Without the stored
        geometry there is nothing to test and every tile is kept.
        """
        if self.center_ra is None or self.radius_deg is None:
            return True
        return _angular_separation(
            ra_deg, dec_deg, self.center_ra, self.center_dec
        ) <= (self.radius_deg + float(radius_deg))

    def __repr__(self):  # pragma: no cover - debugging aid
        return f"IndexEntry({self.name}, {self.scale_lo_arcmin:.1f}-{self.scale_hi_arcmin:.1f}')"


class IndexManifest:
    """The catalogue of available index files."""

    def __init__(self, entries):
        self.entries = list(entries)

    @classmethod
    def load(cls, path=None):
        path = path or MANIFEST_PATH
        with open(path) as handle:
            data = json.load(handle)
        return cls(IndexEntry.from_dict(e) for e in data["entries"])

    @classmethod
    def default(cls):
        """The shipped manifest, or an empty one when it has not been built."""
        try:
            return cls.load()
        except (OSError, ValueError):
            return cls([])

    def save(self, path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w") as handle:
            json.dump(
                {"entries": [e.to_dict() for e in self.entries]}, handle, indent=1
            )

    def __len__(self):
        return len(self.entries)


# --- HEALPix ------------------------------------------------------------------

def _healpix_of(ra_deg, dec_deg, nside):
    """HEALPix id in Astrometry.net's own scheme, or None if unavailable.

    ``astrometry.util.util`` ships with the Astrometry.net suite. It is not
    guaranteed to be installed -- a pip-only machine has the in-process solver
    instead, which does not provide it -- so every caller has a geometric
    fallback.
    """
    try:
        from astrometry.util.util import radecdegtohealpix
    except ImportError:
        return None
    try:
        return int(radecdegtohealpix(float(ra_deg), float(dec_deg), int(nside)))
    except Exception:  # pragma: no cover - depends on the bindings
        return None


def healpix_available():
    """True when the native HEALPix bindings can be used for exact selection."""
    return _healpix_of(0.0, 0.0, 1) is not None


def _cone_samples(ra_deg, dec_deg, radius_deg, n_boundary=720):
    """Points on the boundary of a search cone, plus its centre.

    Sampling the boundary rather than only the centre is what makes a pointing
    near a tile edge pull in its neighbours instead of silently missing them.
    """
    ra = np.radians(float(ra_deg))
    dec = np.radians(float(dec_deg))
    radius = np.radians(float(radius_deg))

    samples = [(float(ra_deg), float(dec_deg))]
    if radius <= 0:
        return samples

    angles = np.linspace(0.0, 2.0 * np.pi, int(n_boundary), endpoint=False)
    sin_r, cos_r = np.sin(radius), np.cos(radius)
    sin_d, cos_d = np.sin(dec), np.cos(dec)

    sin_dec2 = sin_d * cos_r + cos_d * sin_r * np.cos(angles)
    sin_dec2 = np.clip(sin_dec2, -1.0, 1.0)
    dec2 = np.arcsin(sin_dec2)
    d_ra = np.arctan2(
        np.sin(angles) * sin_r * cos_d,
        cos_r - sin_d * np.sin(dec2),
    )
    ra2 = np.degrees(ra + d_ra) % 360.0
    samples.extend(zip(ra2.tolist(), np.degrees(dec2).tolist(), strict=True))
    return samples


def _healpixes_in_cone(ra_deg, dec_deg, radius_deg, nside):
    """The set of HEALPix ids a cone touches, or None when unavailable."""
    if not healpix_available():
        return None
    ids = set()
    for sample_ra, sample_dec in _cone_samples(ra_deg, dec_deg, radius_deg):
        pixel = _healpix_of(sample_ra, sample_dec, nside)
        if pixel is not None:
            ids.add(pixel)
    return ids or None


# --- Selection ----------------------------------------------------------------

def field_size_arcmin(pixel_scale_arcsec, naxis1, naxis2):
    """Field size from the **short** axis, in arcmin.

    Deliberately conservative: using the short axis means the quad-scale window
    never asks for quads larger than the field can actually contain.
    """
    if not pixel_scale_arcsec or pixel_scale_arcsec <= 0:
        return None
    short = min(int(naxis1 or 0), int(naxis2 or 0))
    if short <= 0:
        return None
    return short * float(pixel_scale_arcsec) / 60.0


def required_indexes(ra_deg, dec_deg, pixel_scale_arcsec, naxis1, naxis2,
                     manifest=None, radius_deg=3.0, scale_lo_frac=0.10,
                     scale_hi_frac=1.00, logger=None):
    """The index entries a frame needs, smallest useful set first.

    Returns ``[]`` when the pointing or scale is unknown, which the caller
    should treat as "fall back to whatever local set exists" rather than as an
    error -- a frame with no pointing is a normal thing to encounter.
    """
    logger = logger or get_logger("cassa_integrate")
    manifest = manifest if manifest is not None else IndexManifest.default()
    if not len(manifest):
        return []

    size = field_size_arcmin(pixel_scale_arcsec, naxis1, naxis2)
    if size is None or ra_deg is None or dec_deg is None:
        return []

    lo, hi = size * float(scale_lo_frac), size * float(scale_hi_frac)

    # HEALPix ids per nside used by the manifest, computed once.
    cones = {}
    for entry in manifest.entries:
        if entry.is_all_sky or entry.nside in cones:
            continue
        cones[entry.nside] = _healpixes_in_cone(ra_deg, dec_deg, radius_deg, entry.nside)

    chosen = []
    for entry in manifest.entries:
        # Scale overlap: keep any index whose quad range meets the window.
        if entry.scale_hi_arcmin < lo or entry.scale_lo_arcmin > hi:
            continue
        if entry.is_all_sky:
            chosen.append(entry)
            continue
        ids = cones.get(entry.nside)
        if ids is None:
            # No native HEALPix numbering available (a pip-only install has no
            # astrometry.net bindings, and no other HEALPix library reproduces
            # their numbering). Fall back to the tile geometry stored in the
            # manifest, which is a tight superset -- a handful of neighbouring
            # tiles rather than all 48.
            if entry.reaches(ra_deg, dec_deg, radius_deg):
                chosen.append(entry)
            continue
        if entry.healpix in ids:
            chosen.append(entry)

    # Finest quads first: they are the ones a small field can actually use.
    chosen.sort(key=lambda e: (e.scale_lo_arcmin, e.name))
    total_mb = sum(e.bytes for e in chosen) / 1e6
    logger.debug(
        "Index selection: field %.1f' -> quads %.2f'-%.2f' -> %d file(s), %.0f MB",
        size, lo, hi, len(chosen), total_mb,
    )
    return chosen


# --- Stores -------------------------------------------------------------------

class LocalIndexStore:
    """Index files already on disk. Never downloads anything."""

    def __init__(self, directory):
        self.directory = directory

    def path_of(self, entry):
        candidate = os.path.join(self.directory, entry.name)
        return candidate if os.path.exists(candidate) else None

    def ensure(self, entries, logger=None):
        return [p for p in (self.path_of(e) for e in entries) if p]


class HttpIndexStore:
    """Index files fetched on demand and cached.

    Downloads to ``<name>.part``, verifies the checksum when the manifest has
    one, and renames atomically, so an interrupted run never leaves a truncated
    file that a later solve would silently fail on.
    """

    def __init__(self, base_url, cache_dir, token=None, cache_gb=20.0,
                 local_dirs=(), download=True):
        self.base_url = base_url.rstrip("/") + "/"
        self.cache_dir = cache_dir
        self.token = token
        self.cache_gb = float(cache_gb)
        #: Directories searched before the cache, so an existing local set wins.
        self.local_dirs = [d for d in local_dirs if d and os.path.isdir(d)]
        self.download = bool(download)

    # -- lookup ---------------------------------------------------------------
    def path_of(self, entry):
        """An existing copy of this index, from a local set or the cache."""
        for directory in self.local_dirs:
            candidate = os.path.join(directory, entry.name)
            if os.path.exists(candidate):
                return candidate
        cached = os.path.join(self.cache_dir, entry.name)
        return cached if os.path.exists(cached) else None

    def ensure(self, entries, logger=None):
        """Return local paths for ``entries``, downloading what is missing."""
        logger = logger or get_logger("cassa_integrate")
        paths, missing = [], []
        for entry in entries:
            existing = self.path_of(entry)
            if existing:
                paths.append(existing)
            else:
                missing.append(entry)

        if missing and not self.download:
            logger.warning(
                "%d index file(s) are missing and downloads are disabled. "
                "Fetch them with:  cassa-index-fetch --from-headers <raw dir>",
                len(missing),
            )
            return paths

        if missing:
            total_mb = sum(e.bytes for e in missing) / 1e6
            logger.info(
                "Fetching %d astrometry index file(s) (%.0f MB) into %s",
                len(missing), total_mb, self.cache_dir,
            )
            self._evict_for(sum(e.bytes for e in missing), keep=entries, logger=logger)

        for entry in missing:
            try:
                paths.append(self._download(entry, logger))
            except Exception as exc:
                logger.warning("Could not fetch %s: %s", entry.name, exc)
        return paths

    # -- transfer -------------------------------------------------------------
    def _request(self, url, offset=0):
        request = urllib.request.Request(url)
        if self.token:
            request.add_header("Authorization", f"Bearer {self.token}")
        if offset:
            request.add_header("Range", f"bytes={offset}-")
        return request

    def _download(self, entry, logger):
        os.makedirs(self.cache_dir, exist_ok=True)
        destination = os.path.join(self.cache_dir, entry.name)
        partial = destination + ".part"
        url = self.base_url + entry.url_path

        last_error = None
        for attempt in range(MAX_RETRIES):
            offset = os.path.getsize(partial) if os.path.exists(partial) else 0
            try:
                with urllib.request.urlopen(self._request(url, offset), timeout=60) as response:
                    # A server that ignores Range restarts the file.
                    mode = "ab" if (offset and response.status == 206) else "wb"
                    if mode == "wb":
                        offset = 0
                    with open(partial, mode) as handle:
                        while True:
                            chunk = response.read(CHUNK_BYTES)
                            if not chunk:
                                break
                            handle.write(chunk)
                break
            except (urllib.error.URLError, OSError, TimeoutError) as exc:
                last_error = exc
                if attempt == MAX_RETRIES - 1:
                    raise
                time.sleep(2 ** attempt)
        else:  # pragma: no cover - loop always breaks or raises
            raise last_error

        if entry.sha256 and _sha256(partial) != entry.sha256:
            os.remove(partial)
            raise ValueError(f"checksum mismatch for {entry.name}")

        os.replace(partial, destination)
        logger.info("  fetched %s (%.0f MB)", entry.name, os.path.getsize(destination) / 1e6)
        return destination

    # -- housekeeping ---------------------------------------------------------
    def _evict_for(self, incoming_bytes, keep=(), logger=None):
        """Make room in the cache, never evicting a file this run needs."""
        logger = logger or get_logger("cassa_integrate")
        if not os.path.isdir(self.cache_dir) or self.cache_gb <= 0:
            return
        protected = {e.name for e in keep}
        files = []
        for name in os.listdir(self.cache_dir):
            if name in protected or not name.endswith(".fits"):
                continue
            path = os.path.join(self.cache_dir, name)
            try:
                stat = os.stat(path)
            except OSError:
                continue
            files.append((stat.st_atime, stat.st_size, path))

        limit = self.cache_gb * 1e9
        used = sum(size for _, size, _ in files) + incoming_bytes
        files.sort()  # least recently used first
        for _, size, path in files:
            if used <= limit:
                break
            try:
                os.remove(path)
                used -= size
                logger.info("  evicted %s from the index cache", os.path.basename(path))
            except OSError:
                pass


def _angular_separation(ra1, dec1, ra2, dec2):
    """Great-circle separation between two sky positions, in degrees."""
    lon1, lat1, lon2, lat2 = (np.radians(float(v)) for v in (ra1, dec1, ra2, dec2))
    cos_sep = (np.sin(lat1) * np.sin(lat2)
               + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2))
    return float(np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0))))


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_store(config, logger=None):
    """The index store this configuration implies.

    A populated ``astrometry_index_dir`` is searched first, so an existing
    installation keeps its behaviour and downloads nothing.
    """
    local = config.resolve_astrometry_index_dir()
    return HttpIndexStore(
        base_url=config.phase2.index_url,
        cache_dir=config.resolve_index_cache_dir(),
        token=config.resolve_index_token(),
        cache_gb=config.phase2.index_cache_gb,
        local_dirs=[local],
        download=config.phase2.index_download,
    )
