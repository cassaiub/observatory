"""``cassa-index-fetch``: get the astrometry indexes a field needs, in advance.

The pipeline fetches index files as it solves, so this command is not required.
It exists for the case that motivated the whole feature: preparing a laptop
before taking it somewhere without a network, and knowing what that will cost
before starting the download.

    cassa-index-fetch --ra 339.27 --dec 34.42 --scale 0.598 --size 2048
    cassa-index-fetch --from-headers raw/ --dry-run
"""

from cassa_photometry.astrometry_index import (
    IndexManifest,
    build_store,
    required_indexes,
)
from cassa_photometry.config import load_config
from cassa_photometry.instruments import get_profile
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import find_raw_frames


def _pointings_from_headers(directory, instrument, logger):
    """Pointing, scale and size for every readable frame in a directory."""
    from astropy.io import fits

    from cassa_photometry.phase2_integration.wcs import _as_degrees

    paths = find_raw_frames(directory)
    if not paths:
        logger.warning("No FITS frames found in %s or below it", directory)
        return []

    pointings, skipped = [], 0
    for path in paths:
        try:
            header = fits.getheader(path)
        except Exception:
            skipped += 1
            continue
        ra = _as_degrees(header.get("OBJCTRA"), header.get("RA"), is_ra=True)
        dec = _as_degrees(header.get("OBJCTDEC"), header.get("DEC"), is_ra=False)
        if ra is None or dec is None:
            skipped += 1
            continue
        pointings.append({
            "ra": ra, "dec": dec,
            "scale": instrument.get_pixel_scale(header),
            "naxis1": header.get("NAXIS1"), "naxis2": header.get("NAXIS2"),
            "path": path,
        })

    logger.info("Read %d frame(s) with a pointing; %d had none.", len(pointings), skipped)
    return pointings


def run(ra=None, dec=None, scale=None, size=None, from_headers=None,
        dry_run=False, wide=False, config=None, logger=None):
    """Select and optionally fetch index files. Returns a process exit code."""
    config = config or load_config()
    logger = logger or get_logger("cassa_index_fetch")
    instrument = get_profile(config.instrument, config=config)

    manifest = IndexManifest.default()
    if not len(manifest):
        logger.error(
            "No index manifest is available. Rebuild it with "
            "tools/build_index_manifest.py --from-server."
        )
        return 1

    phase2 = config.phase2
    scale_lo_frac = (phase2.index_scale_lo_frac_wide if wide
                     else phase2.index_scale_lo_frac)

    if from_headers:
        pointings = _pointings_from_headers(from_headers, instrument, logger)
        if not pointings:
            return 1
    else:
        if not scale or not size:
            logger.error("--ra/--dec also need --scale (arcsec/pixel) and --size (pixels).")
            return 1
        pointings = [{"ra": ra, "dec": dec, "scale": scale,
                      "naxis1": size, "naxis2": size, "path": "(command line)"}]

    # Union across every frame: a night can cover more than one field.
    selected, no_scale = {}, 0
    for pointing in pointings:
        if not pointing["scale"]:
            no_scale += 1
            continue
        for entry in required_indexes(
            pointing["ra"], pointing["dec"], pointing["scale"],
            pointing["naxis1"], pointing["naxis2"], manifest=manifest,
            radius_deg=phase2.index_search_radius_deg,
            scale_lo_frac=scale_lo_frac, scale_hi_frac=phase2.index_scale_hi_frac,
            logger=logger,
        ):
            selected[entry.name] = entry

    if no_scale:
        logger.warning(
            "%d frame(s) had no usable plate scale and were skipped. Set "
            "detector.pixel_scale_arcsec in your config, or write SECPIX / "
            "XPIXSZ+FOCALLEN at acquisition.", no_scale,
        )
    if not selected:
        logger.error("Nothing selected: no frame supplied both a pointing and a scale.")
        return 1

    entries = sorted(selected.values(), key=lambda e: (e.scale_lo_arcmin, e.name))
    total_mb = sum(e.bytes for e in entries) / 1e6

    store = build_store(config)
    have = [e for e in entries if store.path_of(e)]
    need = [e for e in entries if not store.path_of(e)]

    print(f"\n{len(entries)} index file(s) needed, {total_mb:,.0f} MB total")
    print(f"  already present: {len(have)}")
    print(f"  to download:     {len(need)}  ({sum(e.bytes for e in need) / 1e6:,.0f} MB)")
    print(f"  cache:           {config.resolve_index_cache_dir()}\n")
    for entry in entries:
        mark = "have" if store.path_of(entry) else " get"
        print(f"  [{mark}] {entry.name:24s} {entry.scale_lo_arcmin:7.1f}-"
              f"{entry.scale_hi_arcmin:7.1f}'  {entry.bytes / 1e6:7.0f} MB")

    if dry_run:
        print("\n--dry-run: nothing downloaded.")
        return 0
    if not need:
        print("\nEverything is already available.")
        return 0

    paths = store.ensure(entries, logger=logger)
    missing = len(entries) - len(paths)
    if missing:
        logger.error("%d file(s) could not be fetched.", missing)
        return 1
    print(f"\nReady: {len(paths)} index file(s) available locally.")
    return 0
