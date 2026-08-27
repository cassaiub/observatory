"""FITS frame I/O, timing, and target-altitude resolution for the monitor."""

import os
import time as _time
from dataclasses import dataclass

import numpy as np
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


@dataclass
class FrameMeta:
    """Per-frame metadata extracted from the FITS header."""

    time: object = None            # astropy.time.Time or None
    exptime: float = None
    coord: object = None           # astropy SkyCoord or None
    altitude_deg: float = None     # if the header states it directly


def read_frame(path):
    """Read a single 2D frame (float32) and its header from a FITS file.

    If the file is a 3D cube with one plane it is squeezed; a multi-plane cube
    raises (use :func:`read_cube` for those).
    """
    with fits.open(path, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        data = np.asarray(data, dtype=np.float32)
    if data.ndim == 3:
        if data.shape[0] == 1:
            data = data[0]
        else:
            raise ValueError(f"{path} is a multi-plane cube; use read_cube().")
    return data, header


def read_cube(path):
    """Read a 3D cube ``(frames, y, x)`` and its header."""
    with fits.open(path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32)
        header = hdul[0].header
    if data.ndim == 2:
        data = data[np.newaxis, ...]
    return data, header


def extract_meta(header):
    """Pull timing, exposure, pointing, and any stated altitude from a header."""
    meta = FrameMeta()

    for key in ("DATE-OBS", "DATE_OBS", "DATE"):
        if key in header and header[key]:
            try:
                meta.time = Time(str(header[key]), format="isot", scale="utc")
                break
            except Exception:
                continue

    for key in ("EXPTIME", "EXPOSURE"):
        if key in header:
            try:
                meta.exptime = float(header[key])
                break
            except (TypeError, ValueError):
                continue

    meta.coord = _parse_radec(header)

    for key in ("OBJCTALT", "CENTALT", "ALTITUDE", "ALT", "EL"):
        if key in header and header[key] not in (None, ""):
            try:
                meta.altitude_deg = float(header[key])
                break
            except (TypeError, ValueError):
                continue

    return meta


def altitude_deg(meta, site):
    """Return the target altitude, from the header if present else computed.

    Computed from the header RA/Dec + observation time + the configured site.
    Returns None if it cannot be determined.
    """
    if meta.altitude_deg is not None:
        return meta.altitude_deg
    if meta.coord is None or meta.time is None:
        return None
    if site.latitude_deg == 0.0 and site.longitude_deg == 0.0:
        return None  # site not configured
    location = EarthLocation(lat=site.latitude_deg * u.deg,
                             lon=site.longitude_deg * u.deg,
                             height=site.elevation_m * u.m)
    altaz = meta.coord.transform_to(AltAz(obstime=meta.time, location=location))
    return float(altaz.alt.deg)


def file_is_stable(path, delay_s=0.5):
    """Return True if the file finished writing at least ``delay_s`` ago.

    Uses a non-blocking modification-time age check (no per-file sleep), so
    catching up on many pre-existing frames stays fast while a file still being
    written -- with a fresh mtime -- is deferred to a later scan.
    """
    try:
        return os.path.getsize(path) > 0 and (_time.time() - os.path.getmtime(path)) >= delay_s
    except OSError:
        return False


def _parse_radec(header):
    """Parse RA/Dec from common keyword conventions into a SkyCoord, or None."""
    ra = next((header[k] for k in ("OBJCTRA", "RA", "CRVAL1") if k in header and header[k] != ""), None)
    dec = next((header[k] for k in ("OBJCTDEC", "DEC", "CRVAL2") if k in header and header[k] != ""), None)
    if ra is None or dec is None:
        return None
    try:  # numeric degrees
        return SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg)
    except (TypeError, ValueError):
        pass
    try:  # sexagesimal strings (RA in hours, Dec in degrees)
        return SkyCoord(str(ra), str(dec), unit=(u.hourangle, u.deg))
    except Exception:
        return None
