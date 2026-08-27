"""Windowed sub-pixel centroiding and per-doublet differential tracking."""

import numpy as np
from scipy.ndimage import median_filter
from photutils.centroids import centroid_com


def _com_once(frame, x0, y0, box):
    """Single-pass background-subtracted centre of mass in a box around (x0, y0)."""
    half = box // 2
    xi, yi = int(round(x0)), int(round(y0))
    y1, y2 = yi - half, yi + half + 1
    x1, x2 = xi - half, xi + half + 1
    if y1 < 0 or x1 < 0 or y2 > frame.shape[0] or x2 > frame.shape[1]:
        return None
    cut = median_filter(frame[y1:y2, x1:x2], size=3)
    sub = np.maximum(cut - np.median(cut), 0)
    if sub.sum() <= 0:
        return None
    cx, cy = centroid_com(sub)
    if not (np.isfinite(cx) and np.isfinite(cy)):
        return None
    return x1 + cx, y1 + cy, float(sub.sum())


def centroid_in_box(frame, x0, y0, box):
    """Two-pass sub-pixel centroid around ``(x0, y0)``.

    Applies a 3x3 median despike and local-background subtraction (as the
    original DIMM script did), but only within the spot's window so many spots
    can be measured independently. A second pass re-centres the box on the
    first-pass centroid, which removes the inward bias that common (tracking)
    motion would otherwise introduce by clipping the spot's wings. Returns
    ``(x, y, flux)`` in global coordinates, or ``None`` if out of bounds / empty.
    """
    first = _com_once(frame, x0, y0, box)
    if first is None:
        return None
    second = _com_once(frame, first[0], first[1], box)
    return second if second is not None else first


def track_doublet(frames, doublet, config):
    """Measure a doublet's per-frame separation across a window of frames.

    Returns ``(dx, dy, n_valid)`` where ``dx``/``dy`` are the per-frame
    ``spot1 - spot2`` separations (pixels) for the frames in which both spots
    were successfully centroided.
    """
    box = config.detection.centroid_box
    dxs, dys = [], []
    for fr in frames:
        c1 = centroid_in_box(fr, doublet.x1, doublet.y1, box)
        c2 = centroid_in_box(fr, doublet.x2, doublet.y2, box)
        if c1 is None or c2 is None:
            continue
        dxs.append(c1[0] - c2[0])
        dys.append(c1[1] - c2[1])
    return np.array(dxs), np.array(dys), len(dxs)
