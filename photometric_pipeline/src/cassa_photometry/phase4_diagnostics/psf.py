"""PSF / image-quality primitives shared by the stage diagnostics.

Stage-independent measurements: detect stars, fit their profiles to estimate the
FWHM and ellipticity, build a stacked radial profile, and summarise basic image
statistics. Uses tools already relied on elsewhere in the package (DAOStarFinder,
Background2D) plus astropy.modeling for the 2D Gaussian fit.
"""

import warnings

import numpy as np
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.modeling import models, fitting
from photutils.detection import DAOStarFinder
from photutils.background import Background2D, MedianBackground

_FWHM_PER_SIGMA = 2.3548200450309493  # 2*sqrt(2*ln2)


def image_stats(data):
    """Return basic robust statistics of an image."""
    finite = data[np.isfinite(data)]
    mean, median, std = sigma_clipped_stats(finite, sigma=3.0)
    dmax = float(np.max(finite)) if finite.size else np.nan
    return {
        "min": float(np.min(finite)) if finite.size else np.nan,
        "max": dmax,
        "median": float(median),
        "mean": float(mean),
        "std": float(std),
        "p1": float(np.percentile(finite, 1)) if finite.size else np.nan,
        "p99": float(np.percentile(finite, 99)) if finite.size else np.nan,
        # Fraction of pixels near the peak value (a rough "hot/saturated" gauge).
        "hot_frac": float(np.mean(finite >= 0.95 * dmax)) if finite.size else np.nan,
    }


def background_rms(data, box=(64, 64)):
    """Return the median 2D-background RMS (falls back to a robust scalar)."""
    try:
        bkg = Background2D(data, box, filter_size=(3, 3), bkg_estimator=MedianBackground())
        return float(np.nanmedian(bkg.background_rms))
    except Exception:
        _, _, std = sigma_clipped_stats(data, sigma=3.0)
        return float(std)


def detect_stars(data, fwhm_guess=3.5, threshold=5.0):
    """Detect point sources with DAOStarFinder; returns a table or None."""
    _, median, std = sigma_clipped_stats(data, sigma=3.0)
    if not np.isfinite(std) or std <= 0:
        return None
    finder = DAOStarFinder(fwhm=fwhm_guess, threshold=threshold * std)
    return finder(np.nan_to_num(data - median))


def _isolated_bright_stars(sources, shape, cutout, max_stars):
    """Pick isolated, non-edge stars, brightest first (saturation handled at fit time)."""
    ny, nx = shape
    order = np.argsort(sources["flux"])[::-1]
    xs = np.asarray(sources["xcentroid"])
    ys = np.asarray(sources["ycentroid"])

    chosen = []
    for i in order:
        x, y = xs[i], ys[i]
        if x < cutout or x > nx - cutout or y < cutout or y > ny - cutout:
            continue
        # Reject if another detected source sits within 2 cutouts.
        d = np.hypot(xs - x, ys - y)
        d[i] = np.inf
        if np.any(d < 2 * cutout):
            continue
        chosen.append((x, y))
        if len(chosen) >= max_stars:
            break
    return chosen


def _fit_gaussian_fwhm(cutout):
    """Fit a 2D Gaussian + constant to a star cutout -> (fwhm_px, ellipticity)."""
    ny, nx = cutout.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    bg = np.median(cutout)
    peak = float(np.max(cutout) - bg)
    if peak <= 0:
        return None
    g = models.Gaussian2D(amplitude=peak, x_mean=nx / 2, y_mean=ny / 2,
                          x_stddev=2.0, y_stddev=2.0) + models.Const2D(bg)
    fitter = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fit = fitter(g, xx, yy, cutout, maxiter=200)
        except Exception:
            return None
    sx, sy = abs(fit.x_stddev_0.value), abs(fit.y_stddev_0.value)
    if not (np.isfinite(sx) and np.isfinite(sy)) or max(sx, sy) > nx:
        return None
    fwhm = _FWHM_PER_SIGMA * np.sqrt(sx * sy)          # geometric-mean FWHM
    ellip = 1.0 - min(sx, sy) / max(sx, sy)
    return fwhm, ellip


def _radial_profile(cutout, nbins=None):
    """Azimuthally-averaged, peak-normalised radial profile of one cutout."""
    ny, nx = cutout.shape
    cy, cx = (ny - 1) / 2, (nx - 1) / 2
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - cx, yy - cy)
    bg = np.median(cutout)
    prof = cutout - bg
    peak = np.max(prof)
    if peak <= 0:
        return None
    prof = prof / peak
    rmax = int(min(cx, cy))
    nbins = nbins or rmax
    bins = np.linspace(0, rmax, nbins + 1)
    idx = np.digitize(r.ravel(), bins) - 1
    vals = prof.ravel()
    radii, profile = [], []
    for b in range(nbins):
        m = idx == b
        if np.any(m):
            radii.append(0.5 * (bins[b] + bins[b + 1]))
            profile.append(np.mean(vals[m]))
    return np.array(radii), np.array(profile)


def estimate_fwhm(data, fwhm_guess=3.5, threshold=5.0, pixscale=None,
                  max_stars=25, cutout=15):
    """Estimate the image FWHM by fitting isolated bright stars.

    Returns
    -------
    dict with keys: ``fwhm_px``, ``fwhm_arcsec`` (or None), ``ellipticity``,
    ``n_detected``, ``n_used``, ``radii``, ``profile`` (stacked radial profile),
    and ``positions`` (list of (x, y) used).
    """
    result = {"fwhm_px": np.nan, "fwhm_arcsec": None, "ellipticity": np.nan,
              "n_detected": 0, "n_used": 0, "radii": None, "profile": None,
              "positions": []}
    sources = detect_stars(data, fwhm_guess, threshold)
    if sources is None or len(sources) == 0:
        return result
    result["n_detected"] = len(sources)

    stars = _isolated_bright_stars(sources, data.shape, cutout, max_stars)
    half = cutout // 2
    fwhms, ellips, profiles, used = [], [], [], []
    for x, y in stars:
        xi, yi = int(round(x)), int(round(y))
        cut = data[yi - half:yi + half + 1, xi - half:xi + half + 1]
        if cut.shape != (cutout, cutout) or not np.all(np.isfinite(cut)):
            continue
        # Skip truly saturated stars: a flat top has many pixels at the peak value.
        if np.sum(cut == np.max(cut)) > 3:
            continue
        fit = _fit_gaussian_fwhm(cut)
        if fit is None:
            continue
        fwhms.append(fit[0])
        ellips.append(fit[1])
        used.append((x, y))
        rp = _radial_profile(cut)
        if rp is not None:
            profiles.append(rp)

    if not fwhms:
        return result
    fwhms = np.array(fwhms)
    clipped = fwhms[~sigma_clip(fwhms, sigma=3.0, maxiters=3).mask]
    fwhm_px = float(np.median(clipped)) if clipped.size else float(np.median(fwhms))

    result["fwhm_px"] = fwhm_px
    result["fwhm_arcsec"] = None if not pixscale else fwhm_px * float(pixscale)
    result["ellipticity"] = float(np.median(ellips))
    result["n_used"] = len(fwhms)
    result["positions"] = used

    if profiles:
        # Interpolate each profile onto a common radial grid, then average.
        rmax = min(p[0][-1] for p in profiles)
        grid = np.linspace(0, rmax, 20)
        stack = [np.interp(grid, r, p) for r, p in profiles]
        result["radii"] = grid
        result["profile"] = np.mean(stack, axis=0)
    return result
