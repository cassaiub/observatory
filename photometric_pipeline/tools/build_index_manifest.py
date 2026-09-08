#!/usr/bin/env python
"""Build the astrometry index manifest shipped with the package.

The manifest is what lets the pipeline decide which index files a field needs
without downloading any of them: one entry per file, giving its sky tile and the
range of quad scales it holds.

Run this once when the index server's contents change, and commit the result:

    python tools/build_index_manifest.py --from-server \\
        -o src/cassa_photometry/data/index_manifest.json

    # or, to take the scales from index files you already hold:
    python tools/build_index_manifest.py --from-dir astrometry_data --from-server \\
        -o src/cassa_photometry/data/index_manifest.json

Scales come from the files' own ``SCALE_L``/``SCALE_U`` headers whenever a local
copy of that series exists, and otherwise from :data:`SERIES_SCALES` below.
"""

import argparse
import glob
import json
import os
import re
import sys
import urllib.request

import numpy as np

#: Quad-scale range per series, in arcmin. Every value here was read from the
#: SCALE_L/SCALE_U headers of a real index file except the three series this
#: observatory does not hold (4200, 4201, 4207), which continue the same
#: geometric progression and are marked. A local file always overrides this.
SERIES_SCALES = {
    4200: (2.0, 2.8),      # inferred
    4201: (2.8, 4.0),      # inferred
    4202: (4.0, 5.6),
    4203: (5.6, 8.0),
    4204: (8.0, 11.0),
    4205: (11.0, 16.0),
    4206: (16.0, 22.0),
    4207: (22.0, 30.0),    # inferred
    4208: (30.0, 42.0), 4209: (42.0, 60.0), 4210: (60.0, 85.0),
    4211: (85.0, 120.0), 4212: (120.0, 170.0), 4213: (170.0, 240.0),
    4214: (240.0, 340.0), 4215: (340.0, 480.0), 4216: (480.0, 680.0),
    4217: (680.0, 1000.0), 4218: (1000.0, 1400.0), 4219: (1400.0, 2000.0),
    4107: (22.0, 30.0), 4108: (30.0, 42.0), 4109: (42.0, 60.0),
    4110: (60.0, 85.0), 4111: (85.0, 120.0), 4112: (120.0, 170.0),
    4113: (170.0, 240.0), 4114: (240.0, 340.0), 4115: (340.0, 480.0),
    4116: (480.0, 680.0), 4117: (680.0, 1000.0), 4118: (1000.0, 1400.0),
    4119: (1400.0, 2000.0),
}

#: Directories on the index server to enumerate, and nothing else.
SERVER_DIRECTORIES = ("4100/", "4200/")

#: File count per series -> HEALPix nside. An all-sky series is one file.
COUNT_TO_NSIDE = {48: 2, 12: 1, 1: 0}

FILE_PATTERN = re.compile(r"^index-(\d{4})(?:-(\d{2}))?\.fits$")


def parse_name(name):
    """``("index-4203-14.fits") -> (4203, 14)``; healpix is -1 for an all-sky file."""
    match = FILE_PATTERN.match(name)
    if not match:
        return None
    series = int(match.group(1))
    healpix = int(match.group(2)) if match.group(2) is not None else -1
    return series, healpix


def scales_from_local(directory):
    """Quad-scale range per series, read from whatever files are on disk."""
    scales = {}
    if not directory or not os.path.isdir(directory):
        return scales
    from astropy.io import fits

    for path in sorted(glob.glob(os.path.join(directory, "index-*.fits"))):
        parsed = parse_name(os.path.basename(path))
        if not parsed or parsed[0] in scales:
            continue
        try:
            header = fits.getheader(path, 0)
            scales[parsed[0]] = (
                float(np.degrees(header["SCALE_L"]) * 60.0),
                float(np.degrees(header["SCALE_U"]) * 60.0),
            )
        except Exception as exc:  # pragma: no cover - depends on the files
            print(f"  ! could not read {os.path.basename(path)}: {exc}", file=sys.stderr)
    return scales


def entries_from_server(base_url, timeout=60):
    """Enumerate every index file the server offers, with its size."""
    found = []
    for directory in SERVER_DIRECTORIES:
        url = base_url.rstrip("/") + "/" + directory
        try:
            html = urllib.request.urlopen(url, timeout=timeout).read().decode("utf8", "replace")
        except Exception as exc:
            print(f"  ! {url}: {exc}", file=sys.stderr)
            continue
        # Apache's autoindex lists "<a href="name">name</a>  date  size".
        for match in re.finditer(
            r'href="(index-[^"]+\.fits)"[^<]*</a>\s+\S+\s+\S+\s+(\d+)', html
        ):
            name, size = match.group(1), int(match.group(2))
            parsed = parse_name(name)
            if parsed:
                found.append((name, parsed[0], parsed[1], size, directory + name))
        print(f"  {url}: {len([f for f in found if f[4].startswith(directory)])} files")
    return found


def entries_from_dir(directory):
    """Enumerate the index files already on disk."""
    found = []
    for path in sorted(glob.glob(os.path.join(directory, "index-*.fits"))):
        name = os.path.basename(path)
        parsed = parse_name(name)
        if parsed:
            found.append((name, parsed[0], parsed[1], os.path.getsize(path), name))
    return found


def tile_geometry(healpix, nside):
    """Centre and angular radius of a HEALPix tile, in degrees.

    Stored in the manifest so that selection works **without** the
    astrometry.net Python bindings. Those bindings use their own tile
    numbering, which neither ``astropy_healpix`` ordering reproduces (verified:
    17/400 and 190/400 agreement at nside 2), so a pip-only installation cannot
    compute the tile ids itself. With a centre and a radius it does not need to:
    a search cone either reaches the tile or it does not.

    Returns ``None`` when the bindings are unavailable, in which case the
    manifest is built without geometry and the fallback stays coarse.
    """
    if healpix < 0 or nside < 1:
        return None
    try:
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astrometry.util.util import healpix_to_radecdeg
    except ImportError:
        return None
    try:
        centre = healpix_to_radecdeg(int(healpix), int(nside), 0.5, 0.5)
        corners = [healpix_to_radecdeg(int(healpix), int(nside), dx, dy)
                   for dx, dy in ((0, 0), (0, 1), (1, 0), (1, 1))]
    except Exception:
        return None

    middle = SkyCoord(centre[0] * u.deg, centre[1] * u.deg)
    radius = max(
        middle.separation(SkyCoord(x * u.deg, y * u.deg)).deg for x, y in corners
    )
    return float(centre[0]), float(centre[1]), float(radius)


def build(found, scales):
    """Turn raw (name, series, healpix, size, path) rows into manifest entries."""
    counts = {}
    for _, series, _, _, _ in found:
        counts[series] = counts.get(series, 0) + 1

    entries = []
    for name, series, healpix, size, url_path in sorted(found):
        scale = scales.get(series) or SERIES_SCALES.get(series)
        if scale is None:
            print(f"  ! no quad scale known for series {series}; skipping {name}",
                  file=sys.stderr)
            continue
        nside = COUNT_TO_NSIDE.get(counts[series], 0) if healpix >= 0 else 0
        entry = {
            "name": name, "series": series, "healpix": healpix, "nside": nside,
            "scale_lo_arcmin": scale[0], "scale_hi_arcmin": scale[1],
            "bytes": size, "sha256": None, "url_path": url_path,
        }
        geometry = tile_geometry(healpix, nside)
        if geometry:
            entry["center_ra"], entry["center_dec"], entry["radius_deg"] = geometry
        entries.append(entry)
    return entries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--from-server", action="store_true",
                        help="Enumerate the public index server (the usual source).")
    parser.add_argument("--url", default="https://data.astrometry.net/",
                        help="Index server base URL.")
    parser.add_argument("--from-dir", default=None,
                        help="Read quad scales (and, without --from-server, the file "
                             "list) from a local index directory.")
    parser.add_argument("-o", "--output",
                        default="src/cassa_photometry/data/index_manifest.json")
    args = parser.parse_args(argv)

    if not args.from_server and not args.from_dir:
        parser.error("give --from-server, --from-dir, or both.")

    scales = scales_from_local(args.from_dir)
    if scales:
        print(f"Quad scales read from {len(scales)} local series.")

    found = entries_from_server(args.url) if args.from_server else entries_from_dir(args.from_dir)
    entries = build(found, scales)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump({
            "source": args.url if args.from_server else os.path.abspath(args.from_dir),
            "entries": entries,
        }, handle, indent=1)

    total_gb = sum(e["bytes"] for e in entries) / 1e9
    tiled = sum(1 for e in entries if e["healpix"] >= 0)
    print(f"\nWrote {len(entries)} entries ({tiled} tiled, {len(entries) - tiled} all-sky), "
          f"{total_gb:.1f} GB total -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
