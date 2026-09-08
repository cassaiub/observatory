"""Bin frames into observing epochs.

Without this, phase 2's stacking key is object/filter/camera/exposure and
carries no date, so five nights of the same field collapse into a single
master -- the time axis is averaged away before phase 3 ever sees it, and no
light curve is possible.

**Binning is on local noon, not the UT calendar date.** An observing night
crosses UT midnight at most longitudes, so grouping on the ``DATE-OBS`` date
string cuts a single night in two and produces two half-depth masters that
differ in nothing but which side of midnight they fell. Local solar time is what
makes "the night of the 2nd" mean one thing, and ``SITELONG`` is what converts
one into the other.
"""

import numpy as np

#: Config value meaning "do not bin at all", reproducing the historical grouping.
NO_EPOCH = "none"

#: Config value meaning "one bin per observing night".
NIGHTLY = "night"


#: Conditions already warned about, so a per-frame problem is reported once.
_WARNED = set()


def epoch_key(meta, epoch_bin=NIGHTLY, logger=None):
    """The epoch label for one frame, or ``""`` when epochs are disabled.

    Parameters
    ----------
    meta : dict
        Frame metadata carrying ``mjd_obs``/``date_obs`` and ``site_long``.
    epoch_bin : str
        ``"none"``, ``"night"``, or a duration such as ``"6h"``/``"30m"``.
    """
    if not epoch_bin or str(epoch_bin).strip().lower() == NO_EPOCH:
        return ""

    mjd = _mjd_of(meta)
    if mjd is None:
        if logger is not None and "no-time" not in _WARNED:
            _WARNED.add("no-time")
            logger.warning(
                "Frames carry no DATE-OBS or MJD-OBS, so they cannot be binned "
                "into epochs; they will all be stacked together. Write DATE-OBS "
                "at acquisition to enable time-domain work."
            )
        return ""

    setting = str(epoch_bin).strip().lower()
    if setting == NIGHTLY:
        return _night_label(mjd, meta.get("site_long"), logger)

    hours = _duration_hours(setting)
    if hours is None:
        if logger is not None and setting not in _WARNED:
            _WARNED.add(setting)
            logger.warning("Unrecognised epoch_bin %r; falling back to nightly.", epoch_bin)
        return _night_label(mjd, meta.get("site_long"), logger)

    # Fixed-width bins, anchored on the same local-noon boundary so a sub-night
    # bin still nests inside its night rather than straddling two.
    offset = _noon_offset(meta.get("site_long"), logger)
    bins_per_day = 24.0 / hours
    index = int(np.floor((mjd + offset) * bins_per_day))
    return f"E{index:d}"


def _night_label(mjd, site_long, logger=None):
    """Label for the observing night an instant falls in."""
    offset = _noon_offset(site_long, logger)
    return f"N{int(np.floor(mjd + offset)):d}"


def _noon_offset(site_long, logger=None):
    """Days to add to an MJD so that local noon falls on an integer boundary.

    Local mean solar time leads UT by ``longitude/360`` of a day. Shifting by a
    further half day puts the boundary at local noon, in the middle of the day,
    so a night's frames land on one side of it.
    """
    if site_long is None or site_long == "":
        if logger is not None and "no-longitude" not in _WARNED:
            _WARNED.add("no-longitude")
            logger.warning(
                "No SITELONG, so epochs are binned on a UT noon boundary rather "
                "than local noon. At longitudes far from Greenwich this can split "
                "one observing night into two epochs."
            )
        return -0.5
    try:
        longitude = float(site_long)
    except (TypeError, ValueError):
        return -0.5
    return longitude / 360.0 - 0.5


def _mjd_of(meta):
    """MJD of a frame, from ``MJD-OBS`` or by parsing ``DATE-OBS``."""
    mjd = meta.get("mjd_obs")
    if mjd not in (None, ""):
        try:
            value = float(mjd)
            if np.isfinite(value):
                return value
        except (TypeError, ValueError):
            pass

    date_obs = meta.get("date_obs")
    if not date_obs:
        return None
    try:
        from astropy.time import Time

        return float(Time(str(date_obs), scale="utc").mjd)
    except Exception:
        return None


def _duration_hours(text):
    """``"6h"`` -> 6.0, ``"90m"`` -> 1.5, ``"3600s"`` -> 1.0; None if unparseable."""
    text = str(text).strip().lower()
    for suffix, factor in (("h", 1.0), ("m", 1.0 / 60.0), ("s", 1.0 / 3600.0), ("d", 24.0)):
        if text.endswith(suffix):
            try:
                return float(text[:-1]) * factor
            except ValueError:
                return None
    try:
        return float(text)  # bare number means hours
    except ValueError:
        return None
