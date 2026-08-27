"""Source detection and prism-doublet pairing (the multi-star core).

A DIMM mask makes every star appear twice (the wedge prism shifts one image by a
fixed vector). This module detects all spots on a stacked reference frame,
auto-calibrates that prism vector from the data, and pairs spots into doublets.
"""

from dataclasses import dataclass

import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder


@dataclass
class Source:
    x: float
    y: float
    flux: float
    peak: float
    snr: float
    ellipticity: float = 0.0


@dataclass
class Doublet:
    x1: float
    y1: float
    x2: float
    y2: float
    snr: float
    sep: float


def detect_sources(image, config):
    """Detect point sources on a (stacked) reference image via DAOStarFinder.

    Applies wide-field quality cuts: an optional central-field restriction
    (``central_radius_arcmin``) and an ellipticity cut (``qc.max_ellipticity``)
    that rejects coma-smeared off-axis spots.
    """
    det = config.detection
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    if not np.isfinite(std) or std <= 0:
        return []
    finder = DAOStarFinder(fwhm=det.fwhm_pixels, threshold=det.snr_threshold * std)
    bg_sub = np.nan_to_num(image - median)
    tbl = finder(bg_sub)
    if tbl is None or len(tbl) == 0:
        return []

    ny, nx = image.shape
    m = det.edge_margin_pixels
    cx, cy = nx / 2.0, ny / 2.0
    central_r = None
    if det.central_radius_arcmin > 0:
        central_r = det.central_radius_arcmin * 60.0 / config.hardware.plate_scale_arcsec
    max_ellip = config.qc.max_ellipticity

    sources = []
    for row in tbl:
        x, y = float(row["xcentroid"]), float(row["ycentroid"])
        if x < m or x > nx - m or y < m or y > ny - m:
            continue
        if central_r is not None and np.hypot(x - cx, y - cy) > central_r:
            continue
        ellip = _stamp_ellipticity(bg_sub, x, y, det.ellipticity_box)
        if np.isfinite(ellip) and ellip > max_ellip:
            continue  # reject elongated / aberrated spots
        peak = float(row["peak"])
        sources.append(Source(x, y, float(row["flux"]), peak, peak / std,
                              ellip if np.isfinite(ellip) else 0.0))
    sources.sort(key=lambda s: s.flux, reverse=True)

    # Deduplicate: a spot broadened by tracking motion can be split into several
    # detections by DAOStarFinder. Keep the brightest and drop any within
    # 2*FWHM (well below a realistic prism separation of tens of pixels).
    dedup = 2.0 * det.fwhm_pixels
    kept = []
    for s in sources:
        if all(np.hypot(s.x - k.x, s.y - k.y) > dedup for k in kept):
            kept.append(s)
    return kept[: det.max_sources]


def _stamp_ellipticity(image, x, y, box):
    """Ellipticity ``1 - b/a`` from intensity-weighted second moments in a stamp."""
    half = box // 2
    xi, yi = int(round(x)), int(round(y))
    y1, y2, x1, x2 = yi - half, yi + half + 1, xi - half, xi + half + 1
    if y1 < 0 or x1 < 0 or y2 > image.shape[0] or x2 > image.shape[1]:
        return np.nan
    cut = np.clip(image[y1:y2, x1:x2].astype(float) - np.median(image[y1:y2, x1:x2]), 0, None)
    total = cut.sum()
    if total <= 0:
        return np.nan
    yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
    mx = (xx * cut).sum() / total
    my = (yy * cut).sum() / total
    mxx = ((xx - mx) ** 2 * cut).sum() / total
    myy = ((yy - my) ** 2 * cut).sum() / total
    mxy = ((xx - mx) * (yy - my) * cut).sum() / total
    common = 0.5 * (mxx + myy)
    diff = np.sqrt(max(((mxx - myy) / 2.0) ** 2 + mxy ** 2, 0.0))
    l1, l2 = common + diff, common - diff
    if l1 <= 0:
        return np.nan
    return float(1.0 - np.sqrt(max(l2, 0.0)) / np.sqrt(l1))


def estimate_prism_vector(sources, config):
    """Auto-calibrate the doublet separation vector from repeated pairwise offsets.

    Returns ``(vec, n_support)`` where ``vec`` is the (dx, dy) prism shift in
    pixels, or ``None`` if no consistent offset is found. A configured
    ``hardware.prism_dx_pix/prism_dy_pix`` overrides the auto-detection.
    """
    hw = config.hardware
    if hw.prism_dx_pix is not None and hw.prism_dy_pix is not None:
        return np.array([hw.prism_dx_pix, hw.prism_dy_pix], dtype=float), -1

    tol = config.detection.pair_tol_pixels
    pos = np.array([[s.x, s.y] for s in sources])
    n = len(pos)
    if n < 2:
        return None

    offs = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = pos[j] - pos[i]
            if np.hypot(*d) < 3.0:
                continue
            if d[0] < 0 or (d[0] == 0 and d[1] < 0):   # canonical half-plane
                d = -d
            offs.append(d)
    if not offs:
        return None
    offs = np.array(offs)

    best_i, best_count = 0, -1
    for i in range(len(offs)):
        cnt = int(np.sum(np.hypot(offs[:, 0] - offs[i, 0], offs[:, 1] - offs[i, 1]) < tol))
        if cnt > best_count:
            best_count, best_i = cnt, i
    mask = np.hypot(offs[:, 0] - offs[best_i, 0], offs[:, 1] - offs[best_i, 1]) < tol
    vec = offs[mask].mean(axis=0)
    return vec, int(mask.sum())


def pair_doublets(sources, prism_vec, config):
    """Greedily pair spots (brightest first) whose separation matches the prism vector."""
    tol = config.detection.pair_tol_pixels
    pos = np.array([[s.x, s.y] for s in sources])
    used = [False] * len(sources)
    doublets = []

    for i in range(len(sources)):        # sources are pre-sorted brightest first
        if used[i]:
            continue
        best_j, best_dist = -1, tol
        for j in range(len(sources)):
            if j == i or used[j]:
                continue
            d = pos[j] - pos[i]
            dist = min(np.hypot(*(d - prism_vec)), np.hypot(*(d + prism_vec)))
            if dist < best_dist:
                best_dist, best_j = dist, j
        if best_j >= 0:
            used[i] = used[best_j] = True
            snr = 0.5 * (sources[i].snr + sources[best_j].snr)
            sep = float(np.hypot(*(pos[best_j] - pos[i])))
            doublets.append(Doublet(pos[i][0], pos[i][1], pos[best_j][0], pos[best_j][1], snr, sep))
    return doublets
