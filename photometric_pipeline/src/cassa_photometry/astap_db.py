"""On-demand ASTAP star-database tiles: work out what a field needs, fetch only that.

ASTAP ships its star databases as one large ZIP (859 MB for ``d50``) and offers
no per-tile URL. Downloading the whole thing to solve one field is the same
waste the astrometry.net index selector already avoids, and it breaks the moment
a user repoints: there is nothing to top up incrementally.

Two facts make incremental fetching possible anyway:

* the download host honours HTTP **range** requests, and
* a ZIP is random-access by design -- its central directory lists every member's
  offset and compressed size, so a member can be read without the rest.

So this module treats the published archive as a remote filesystem. Measured on
a real field (NGC7331, 10' across): the central directory costs 3 requests and
85 KB, the 10 tiles it needs cost 6.2 MB, and the solve then takes 0.1 s --
against 859 MB for the full archive, a 138x saving, with nothing to host.

**The tiling.** ASTAP divides the sky into 36 declination bands of 5 degrees,
each holding a band-dependent number of RA tiles -- 1 at the poles, 69 at the
equator. A tile is named ``<series>_<BB><TT>.<ext>``: band, then RA index, both
1-based and zero-padded. The first tile doubles as the "database present"
marker: without ``<series>_0101`` ASTAP reports *"no star database found"*
however many other tiles are on disk, which is why it is always fetched.
"""

import io
import os
import urllib.request
import zipfile

#: RA tiles per declination band, from the south pole to the north. Symmetric,
#: and the sum is the tile count of every series (1476 for d50 and d80).
TILES_PER_BAND = (
    1, 3, 9, 15, 21, 27, 33, 38, 43, 48, 52, 56, 60, 63, 65, 67, 68, 69,
    69, 68, 67, 65, 63, 60, 56, 52, 48, 43, 38, 33, 27, 21, 15, 9, 3, 1,
)

#: Degrees of declination per band.
BAND_HEIGHT_DEG = 180.0 / len(TILES_PER_BAND)

#: Series the pipeline can fetch incrementally, smallest field first. Only
#: series published as a ZIP can be range-fetched: ``d80``'s .deb is xz (not
#: seekable per member) and its .exe is a Windows self-extractor, so ``d80`` is
#: usable only when the user has installed it themselves.
SERIES = {
    "d50": {
        "archive": "d50_star_database.zip",
        # Documented as 6 deg > FOV > 0.2 deg, but the lower bound is soft:
        # verified solving a 0.17 deg field, which is narrower than every CASSA
        # telescope's, so this one series covers the whole fleet.
        "fov_max_deg": 6.0,
        "fov_min_deg": 0.10,
    },
    "d05": {
        "archive": "d05_star_database.zip",
        "fov_max_deg": 6.0,
        "fov_min_deg": 0.6,
    },
    "g05": {
        "archive": "g05_star_database.zip",
        "fov_max_deg": 20.0,
        "fov_min_deg": 3.0,
    },
}

DEFAULT_BASE_URL = ("https://sourceforge.net/projects/astap-program/files/"
                    "star_databases/")


class AstapDatabaseError(RuntimeError):
    """The tiles a field needs could not be obtained."""


def band_of(dec_deg):
    """1-based declination band holding this declination."""
    band = int((float(dec_deg) + 90.0) / BAND_HEIGHT_DEG) + 1
    return min(max(band, 1), len(TILES_PER_BAND))


def tile_of(ra_deg, band):
    """1-based RA tile index within a band."""
    count = TILES_PER_BAND[band - 1]
    step = 360.0 / count
    return int((float(ra_deg) % 360.0) / step) + 1


def series_for_fov(fov_deg):
    """The series covering a field of this height, or None.

    Chooses the densest series whose range contains the field, because a
    too-sparse database is what produces "no solution found" on a narrow field.
    """
    if not fov_deg or fov_deg <= 0:
        return "d50"
    for name in ("d50", "d05", "g05"):
        spec = SERIES[name]
        if spec["fov_min_deg"] <= fov_deg <= spec["fov_max_deg"]:
            return name
    return None


def required_tiles(ra_deg, dec_deg, series="d50", neighbours=1):
    """Tile basenames covering a pointing, plus the marker tile.

    ``neighbours`` widens the selection by that many tiles in each direction, so
    a field near a tile edge -- or a pointing that is a few arcminutes off, which
    is normal -- still has coverage. RA wraps; declination clamps at the poles.
    """
    if ra_deg is None or dec_deg is None:
        return []

    wanted = {f"{series}_0101"}
    centre_band = band_of(dec_deg)
    for band in range(centre_band - neighbours, centre_band + neighbours + 1):
        if not 1 <= band <= len(TILES_PER_BAND):
            continue
        count = TILES_PER_BAND[band - 1]
        centre = tile_of(ra_deg, band)
        for offset in range(-neighbours, neighbours + 1):
            index = (centre - 1 + offset) % count + 1
            wanted.add(f"{series}_{band:02d}{index:02d}")
    return sorted(wanted)


class _HttpRangeFile(io.RawIOBase):
    """A seekable read-only file over HTTP range requests."""

    def __init__(self, url, timeout=120):
        self.url = url
        self.timeout = timeout
        self._pos = 0
        self._size = None
        self.requests = 0
        self.bytes_read = 0

    @property
    def size(self):
        if self._size is None:
            request = urllib.request.Request(self.url, method="HEAD")
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                # Keep the mirror the redirect chose; re-resolving per range
                # request would double the round trips.
                self.url = response.geturl()
                length = response.headers.get("Content-Length")
                if not length:
                    raise AstapDatabaseError(
                        f"{self.url}: no Content-Length, cannot range-fetch.")
                if response.headers.get("Accept-Ranges", "").lower() != "bytes":
                    raise AstapDatabaseError(
                        f"{self.url}: host does not accept range requests.")
                self._size = int(length)
        return self._size

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self._pos

    def seek(self, offset, whence=io.SEEK_SET):
        if whence == io.SEEK_SET:
            self._pos = offset
        elif whence == io.SEEK_CUR:
            self._pos += offset
        else:
            self._pos = self.size + offset
        return self._pos

    def read(self, size=-1):
        if size < 0:
            size = self.size - self._pos
        end = min(self._pos + size, self.size) - 1
        if size == 0 or self._pos > end:
            return b""
        request = urllib.request.Request(
            self.url, headers={"Range": f"bytes={self._pos}-{end}"})
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            data = response.read()
        self.requests += 1
        self.bytes_read += len(data)
        self._pos += len(data)
        return data


class AstapTileStore:
    """Local ASTAP tiles, topped up from the published archive on demand."""

    def __init__(self, directory, base_url=DEFAULT_BASE_URL, download=True,
                 timeout=120, logger=None):
        self.directory = directory
        self.base_url = base_url.rstrip("/") + "/"
        self.download = download
        self.timeout = timeout
        self.logger = logger

    # -- local -----------------------------------------------------------------
    def present(self, stems):
        """Which of these tile basenames are already on disk."""
        if not os.path.isdir(self.directory):
            return set()
        have = {}
        for name in os.listdir(self.directory):
            stem = name.rsplit(".", 1)[0]
            have[stem] = name
        return {s for s in stems if s in have}

    def missing(self, stems):
        return [s for s in stems if s not in self.present(stems)]

    # -- remote ----------------------------------------------------------------
    def _archive_url(self, series):
        spec = SERIES.get(series)
        if spec is None:
            raise AstapDatabaseError(
                f"No incrementally fetchable archive for series {series!r}. "
                f"Known: {', '.join(sorted(SERIES))}.")
        return f"{self.base_url}{spec['archive']}/download"

    def ensure(self, ra_deg, dec_deg, series="d50", neighbours=1):
        """Make sure this field's tiles are on disk. Returns the directory.

        Tiles already present are never re-downloaded, so repointing costs only
        the tiles the new field adds.
        """
        stems = required_tiles(ra_deg, dec_deg, series, neighbours)
        if not stems:
            return self.directory

        absent = self.missing(stems)
        if not absent:
            if self.logger:
                self.logger.info(
                    "    -> ASTAP database: %d tile(s) already cached.", len(stems))
            return self.directory

        if not self.download:
            raise AstapDatabaseError(
                f"{len(absent)} ASTAP tile(s) missing from {self.directory} and "
                f"downloading is disabled (phase2.astap_db_download). "
                f"Needed: {', '.join(absent[:6])}"
                f"{' ...' if len(absent) > 6 else ''}")

        os.makedirs(self.directory, exist_ok=True)
        self._fetch(absent, series)
        return self.directory

    def _fetch(self, stems, series):
        url = self._archive_url(series)
        if self.logger:
            self.logger.info(
                "    -> ASTAP database: fetching %d missing tile(s) from the "
                "published archive (range requests, not the whole file).",
                len(stems))

        backing = _HttpRangeFile(url, timeout=self.timeout)
        try:
            archive = zipfile.ZipFile(backing)
        except Exception as exc:
            raise AstapDatabaseError(
                f"Could not read the ASTAP archive at {url}: {exc}. "
                f"Install a star database manually and set "
                f"phase2.astap_db_dir, or set phase2.astap_db_download: false."
            ) from exc

        members = {}
        for name in archive.namelist():
            members.setdefault(os.path.basename(name).rsplit(".", 1)[0], name)

        written = 0
        for stem in stems:
            member = members.get(stem)
            if member is None:
                # A tile the archive does not contain is not an error: the
                # neighbour expansion deliberately over-asks near the poles.
                continue
            data = archive.read(member)
            target = os.path.join(self.directory, os.path.basename(member))
            temporary = target + ".part"
            with open(temporary, "wb") as handle:
                handle.write(data)
            os.replace(temporary, target)
            written += 1

        if self.logger:
            self.logger.info(
                "    -> ASTAP database: %d tile(s), %.1f MB in %d range request(s).",
                written, backing.bytes_read / 1048576.0, backing.requests)
        if written == 0:
            raise AstapDatabaseError(
                f"None of the required tiles were found in {url}.")


def resolve_db_dir(config):
    """Where ASTAP tiles live.

    Config value > ``CASSA_ASTAP_DB`` > ``~/.cache/cassa-photometry/astap``.
    Mirrors how the astrometry.net index cache resolves, so there is one habit
    to learn rather than two.
    """
    configured = getattr(config.phase2, "astap_db_dir", None)
    if configured:
        return os.path.abspath(os.path.expanduser(configured))
    env = os.environ.get("CASSA_ASTAP_DB")
    if env:
        return os.path.abspath(os.path.expanduser(env))
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(
        os.path.expanduser("~"), ".cache")
    return os.path.join(base, "cassa-photometry", "astap")
