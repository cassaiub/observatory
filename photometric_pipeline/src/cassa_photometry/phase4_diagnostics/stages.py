"""Per-stage diagnostics for the pipeline outputs.

Each ``diagnose_*`` reads one stage's product(s), computes quality metrics, and
builds a diagnostic figure page. All are defensive: missing planes, missing WCS
or too-few stars degrade gracefully rather than raising.
"""

import os

import numpy as np
from astropy.wcs import WCS
from astropy.table import Table

from cassa_photometry.fits_utils import (
    read_mef, DQ_SATURATED, DQ_BAD_PIXEL, DQ_COSMIC_RAY, DQ_NO_DATA,
)
from cassa_photometry.phase2_integration.math_utils import MathEngine
from cassa_photometry.phase4_diagnostics import plots
from cassa_photometry.phase4_diagnostics.psf import (
    estimate_fwhm, image_stats, background_rms,
)

_DQ_BITS = [
    ("SAT", DQ_SATURATED), ("BAD", DQ_BAD_PIXEL),
    ("CR", DQ_COSMIC_RAY), ("NODATA", DQ_NO_DATA),
]


def _pixscale_from_header(header):
    val = header.get("PIXSCALE")
    if val:
        try:
            return float(val)
        except (TypeError, ValueError):
            pass
    try:
        w = WCS(header)
        if w.has_celestial:
            return float(np.abs(w.pixel_scale_matrix[0, 0]) * 3600.0)
    except Exception:
        pass
    return None


def _fwhm_line(fw):
    s = f"FWHM        : {fw['fwhm_px']:.2f} px"
    if fw["fwhm_arcsec"]:
        s += f"  ({fw['fwhm_arcsec']:.2f}\")"
    return s


# --------------------------------------------------------------------------- #
def diagnose_raw(path, config, logger):
    """stage_0 — raw frame sanity + PSF."""
    sci, _, _, header = read_mef(path)
    cfg = config.phase3
    pixscale = _pixscale_from_header(header)
    stats = image_stats(sci)
    fw = estimate_fwhm(sci, cfg.fwhm, cfg.detection_threshold, pixscale)
    sub, _ = MathEngine.extract_2d_background(sci, config.phase2.background_box,
                                              config.phase2.background_filter)
    hot = (sci >= 0.95 * stats["max"]).astype(float)

    fig, ax = plots.new_page(2, 3, f"Stage 0 - Raw: {os.path.basename(path)}")
    im = plots.zscale_panel(ax[0], sci, "Raw image + detected stars")
    plots.overlay_positions(ax[0], fw["positions"])
    plots.histogram_panel(ax[1], sci, "Pixel histogram")
    plots.radial_profile_panel(ax[2], fw["radii"], fw["profile"], fw["fwhm_px"], "Stacked PSF profile")
    plots.zscale_panel(ax[3], sub, "Background-subtracted")
    plots.image_panel(ax[4], hot, "Near-saturation map", cmap="Reds", vmin=0, vmax=1)
    plots.text_panel(ax[5], [
        f"OBJECT      : {header.get('OBJECT', '?')}",
        f"FILTER      : {header.get('FILTER', '?')}",
        f"median      : {stats['median']:.1f}",
        f"std (robust): {stats['std']:.1f}",
        f"max         : {stats['max']:.0f}",
        f"near-sat %  : {100*stats['hot_frac']:.3f}",
        f"stars found : {fw['n_detected']} (used {fw['n_used']})",
        _fwhm_line(fw),
        f"ellipticity : {fw['ellipticity']:.3f}",
    ], "Metrics")

    metrics = {"median": stats["median"], "std": stats["std"], "max": stats["max"],
               "near_sat_frac": stats["hot_frac"], "n_stars": fw["n_used"],
               "fwhm_px": fw["fwhm_px"], "fwhm_arcsec": fw["fwhm_arcsec"],
               "ellipticity": fw["ellipticity"]}
    return metrics, [fig]


# --------------------------------------------------------------------------- #
def diagnose_calibrated(path, config, logger, raw_path=None):
    """stage_1 — calibrated SCI/ERR/DQ, error-plane and DQ sanity."""
    sci, err, dq, header = read_mef(path)
    cfg = config.phase3
    pixscale = _pixscale_from_header(header)
    stats = image_stats(sci)
    fw = estimate_fwhm(sci, cfg.fwhm, cfg.detection_threshold, pixscale)

    dq_counts = {name: int(np.sum((dq & bit) > 0)) for name, bit in _DQ_BITS} if dq is not None else {}
    err_med = float(np.nanmedian(err)) if err is not None else np.nan

    # Background flatness improvement vs the matching raw frame.
    flat_txt = "raw not provided"
    bg_improve = None
    cal_bg = background_rms(sci, config.phase2.background_box)
    if raw_path and os.path.exists(raw_path):
        raw_sci, _, _, _ = read_mef(raw_path)
        raw_bg = background_rms(raw_sci, config.phase2.background_box)
        bg_improve = float(raw_bg / cal_bg) if cal_bg > 0 else None
        flat_txt = f"raw {raw_bg:.1f} -> cal {cal_bg:.1f}"

    fig, ax = plots.new_page(2, 3, f"Stage 1 - Calibrated: {os.path.basename(path)}")
    plots.zscale_panel(ax[0], sci, "SCI [electrons]")
    if err is not None:
        plots.image_panel(ax[1], err, "ERR [electrons]", cmap="magma",
                          vmin=np.nanpercentile(err, 1), vmax=np.nanpercentile(err, 99))
    else:
        plots._text_only(ax[1], "No ERR plane", "ERR")
    if dq is not None:
        plots.image_panel(ax[2], (dq > 0).astype(float), "DQ flagged pixels", cmap="Reds", vmin=0, vmax=1)
        plots.bar_panel(ax[3], [n for n, _ in _DQ_BITS], [dq_counts[n] for n, _ in _DQ_BITS], "DQ flag counts")
    else:
        plots._text_only(ax[2], "No DQ plane", "DQ")
        plots._text_only(ax[3], "No DQ plane", "DQ counts")
    # ERR vs signal (Poisson check)
    if err is not None:
        _err_vs_signal(ax[4], sci, err)
    else:
        plots.radial_profile_panel(ax[4], fw["radii"], fw["profile"], fw["fwhm_px"])
    plots.text_panel(ax[5], [
        f"FILTER      : {header.get('FILTER', '?')}",
        f"BUNIT       : {header.get('BUNIT', '?')}",
        f"median SCI  : {stats['median']:.1f}",
        f"median ERR  : {err_med:.2f}",
        f"bg flatten  : {flat_txt}",
        (f"bg improve  : {bg_improve:.2f}x" if bg_improve else "bg improve  : n/a"),
        f"DQ SAT/BAD  : {dq_counts.get('SAT','-')}/{dq_counts.get('BAD','-')}",
        f"DQ CR/NODATA: {dq_counts.get('CR','-')}/{dq_counts.get('NODATA','-')}",
        _fwhm_line(fw),
    ], "Metrics")

    metrics = {"median_sci": stats["median"], "median_err": err_med,
               "dq_counts": dq_counts, "bg_rms": cal_bg, "bg_improvement": bg_improve,
               "fwhm_px": fw["fwhm_px"], "n_stars": fw["n_used"]}
    return metrics, [fig]


# --------------------------------------------------------------------------- #
def diagnose_master(path, config, logger, single_calibrated_path=None):
    """stage_2 — WCS-solved master stack: depth boost, PSF, WCS status."""
    sci, err, dq, header = read_mef(path)
    cfg = config.phase3
    pixscale = _pixscale_from_header(header)
    stats = image_stats(sci)
    fw = estimate_fwhm(sci, cfg.fwhm, cfg.detection_threshold, pixscale)

    try:
        wcs = WCS(header)
        has_wcs = bool(wcs.has_celestial)
    except Exception:
        has_wcs = False
    stackcnt = header.get("STACKCNT")
    master_noise = background_rms(sci, config.phase2.background_box)
    # Phase 2 records the astrometric fit residual; absent on older masters.
    astrom_rms = header.get("ASTRMS")
    astrom_nstars = header.get("ASTNSTAR")

    snr_boost, expected = None, None
    if single_calibrated_path and os.path.exists(single_calibrated_path):
        s_sci, _, _, _ = read_mef(single_calibrated_path)
        single_noise = background_rms(s_sci, config.phase2.background_box)
        snr_boost = float(single_noise / master_noise) if master_noise > 0 else None
    if stackcnt:
        expected = float(np.sqrt(float(stackcnt)))

    center = "n/a"
    if has_wcs:
        ny, nx = sci.shape
        c = wcs.pixel_to_world(nx / 2, ny / 2)
        center = f"{c.ra.deg:.4f}, {c.dec.deg:.4f}"

    snr_map = np.where((err is not None) & np.isfinite(err) & (err > 0), sci / err, np.nan) \
        if err is not None else None

    fig, ax = plots.new_page(2, 3, f"Stage 2 - Master stack: {os.path.basename(path)}")
    plots.zscale_panel(ax[0], sci, f"Master SCI ({stackcnt or '?'} frames)")
    if err is not None:
        plots.image_panel(ax[1], err, "ERR [electrons]", cmap="magma",
                          vmin=np.nanpercentile(err, 1), vmax=np.nanpercentile(err, 99))
    else:
        plots._text_only(ax[1], "No ERR plane", "ERR")
    if snr_map is not None:
        plots.image_panel(ax[2], snr_map, "SNR map (SCI/ERR)", cmap="viridis",
                          vmin=0, vmax=np.nanpercentile(snr_map, 99))
    else:
        plots._text_only(ax[2], "No ERR plane", "SNR map")
    plots.radial_profile_panel(ax[3], fw["radii"], fw["profile"], fw["fwhm_px"], "Stacked PSF profile")
    if dq is not None:
        plots.image_panel(ax[4], (dq > 0).astype(float), "DQ flagged", cmap="Reds", vmin=0, vmax=1)
    else:
        plots._text_only(ax[4], "No DQ plane", "DQ")
    plots.text_panel(ax[5], [
        f"WCS solved  : {'YES' if has_wcs else 'NO'}",
        (f"astrom RMS  : {astrom_rms:.3f} \"  ({astrom_nstars} stars)"
         if astrom_rms is not None else "astrom RMS  : n/a"),
        f"field centre: {center}",
        f"pixscale    : {pixscale:.3f} \"/px" if pixscale else "pixscale    : n/a",
        f"stack count : {stackcnt or '?'}",
        f"TOT_EXP     : {header.get('TOT_EXP', '?')} s",
        f"master noise: {master_noise:.2f}",
        (f"SNR boost   : {snr_boost:.2f}x (exp ~{expected:.2f}x)"
         if snr_boost and expected else
         (f"SNR boost   : {snr_boost:.2f}x" if snr_boost else "SNR boost   : n/a")),
        _fwhm_line(fw),
    ], "Metrics")

    metrics = {"wcs_solved": has_wcs, "field_center": center, "pixscale": pixscale,
               "astrometric_rms_arcsec": astrom_rms, "astrometry_nstars": astrom_nstars,
               "stack_count": stackcnt, "master_noise": master_noise,
               "snr_boost": snr_boost, "expected_boost": expected,
               "fwhm_px": fw["fwhm_px"], "fwhm_arcsec": fw["fwhm_arcsec"]}
    return metrics, [fig]


# --------------------------------------------------------------------------- #
def diagnose_photometry(catalog_csv, config, logger, fluxcal_fits=None):
    """stage_3 — catalog: photometric errors, depth, ZP, star/galaxy split."""
    cat = Table.read(catalog_csv, format="csv")
    cfg = config.phase3

    zp = zperr = nzp = None
    if fluxcal_fits and os.path.exists(fluxcal_fits):
        _, _, _, header = read_mef(fluxcal_fits)
        zp = header.get("MAGZERO")
        zperr = header.get("MAGZERR")
        nzp = header.get("NZPSTARS")

    mag = np.asarray(cat["Absolute_Mag"], float)
    magerr = np.asarray(cat["Mag_Error"], float) if "Mag_Error" in cat.colnames else np.full_like(mag, np.nan)
    snr = np.asarray(cat["SNR"], float) if "SNR" in cat.colnames else np.full_like(mag, np.nan)
    good = np.isfinite(mag)

    lim_mag = _limiting_mag(mag[good], snr[good], target=5.0)
    turnover = _histogram_turnover(mag[good])

    fig, ax = plots.new_page(2, 3, f"Stage 3 - Photometry: {os.path.basename(catalog_csv)}")
    plots.scatter_panel(ax[0], mag[good], magerr[good], "Magnitude", "Mag error", "Photometric error vs mag")
    plots.scatter_panel(ax[1], mag[good], snr[good], "Magnitude", "SNR", "SNR vs mag", logy=True)
    ax[1].axhline(5, color="red", ls="--", lw=0.8)
    ax[2].hist(mag[good], bins=30, color="steelblue")
    if turnover:
        ax[2].axvline(turnover, color="red", ls="--", lw=0.8, label=f"turnover {turnover:.1f}")
        ax[2].legend(fontsize=8)
    ax[2].set_title("Number counts", fontsize=10)
    ax[2].set_xlabel("Magnitude")
    ax[2].set_ylabel("N")
    if "X_pix" in cat.colnames:
        sc = plots.scatter_panel(ax[3], cat["X_pix"][good], cat["Y_pix"][good],
                                 "X (px)", "Y (px)", "Source map (colour=mag)", c=mag[good])
        fig.colorbar(sc, ax=ax[3], fraction=0.046)
    if "Ellipticity" in cat.colnames:
        plots.scatter_panel(ax[4], mag[good], np.asarray(cat["Ellipticity"], float)[good],
                            "Magnitude", "Ellipticity", "Shape vs mag")
        ax[4].axhline(cfg.ellipticity_star_max, color="green", ls="--", lw=0.8)
    plots.text_panel(ax[5], [
        f"N sources   : {int(np.sum(good))}",
        f"ZP (MAGZERO): {zp:.4f}" if zp is not None else "ZP          : n/a",
        f"ZP err      : {zperr:.4f}" if zperr is not None else "ZP err      : n/a",
        f"N ZP stars  : {nzp}" if nzp is not None else "N ZP stars  : n/a",
        f"limiting mag: {lim_mag:.2f} (SNR=5)" if lim_mag else "limiting mag: n/a",
        f"turnover mag: {turnover:.2f}" if turnover else "turnover mag: n/a",
        f"median magerr: {np.nanmedian(magerr[good]):.3f}",
    ], "Metrics")

    metrics = {"n_sources": int(np.sum(good)), "zero_point": _f(zp), "zero_point_err": _f(zperr),
               "n_zp_stars": _f(nzp), "limiting_mag": lim_mag, "turnover_mag": turnover,
               "median_mag_err": float(np.nanmedian(magerr[good])) if np.any(good) else None}
    return metrics, [fig]


# --------------------------------------------------------------------------- #
def _err_vs_signal(ax, sci, err, nsample=5000):
    """Scatter sampled ERR vs SCI with the sqrt(signal) Poisson expectation."""
    s = sci.ravel()
    e = err.ravel()
    m = np.isfinite(s) & np.isfinite(e) & (s > 0)
    if np.sum(m) > nsample:
        idx = np.random.choice(np.where(m)[0], nsample, replace=False)
    else:
        idx = np.where(m)[0]
    ax.scatter(s[idx], e[idx], s=3, alpha=0.3, color="steelblue")
    xs = np.linspace(0, np.nanpercentile(s[m], 99), 100)
    ax.plot(xs, np.sqrt(np.clip(xs, 0, None)), "r--", lw=1, label="sqrt(signal)")
    ax.set_xlabel("SCI (electrons)")
    ax.set_ylabel("ERR (electrons)")
    ax.set_title("Error vs signal", fontsize=10)
    ax.legend(fontsize=8)


def _limiting_mag(mag, snr, target=5.0, binsize=0.5):
    """Faintest magnitude bin whose median SNR still exceeds ``target``."""
    m = np.isfinite(mag) & np.isfinite(snr) & (snr > 0)
    if np.sum(m) < 5:
        return None
    mag, snr = mag[m], snr[m]
    bins = np.arange(np.floor(mag.min()), np.ceil(mag.max()) + binsize, binsize)
    last = None
    for lo in bins[:-1]:
        sel = (mag >= lo) & (mag < lo + binsize)
        if np.sum(sel) >= 3:
            if np.median(snr[sel]) >= target:
                last = lo + binsize / 2
            else:
                break
    return float(last) if last is not None else None


def _histogram_turnover(mag, binsize=0.5):
    """Magnitude of the number-count peak (a rough completeness limit)."""
    m = np.isfinite(mag)
    if np.sum(m) < 5:
        return None
    counts, edges = np.histogram(mag[m], bins=30)
    return float(0.5 * (edges[np.argmax(counts)] + edges[np.argmax(counts) + 1]))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
