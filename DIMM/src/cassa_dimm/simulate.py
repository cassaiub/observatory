"""Physically-calibrated DIMM simulator.

Given a target seeing, it computes the *expected* differential image-motion
variance from the same Sarazin & Roddier relation the reducer inverts, then
injects it into ``n_stars`` prism doublets. A run therefore has a known answer,
so the pipeline can be validated (recovered seeing ~= injected seeing).
"""

import os
from datetime import datetime, timedelta

import numpy as np
from astropy.io import fits

from cassa_dimm.config import load_config
from cassa_dimm.seeing import response_coefficients, _ARCSEC_TO_RAD, _SEEING_CONST


def _render_gaussian(grid_y, grid_x, x0, y0, sigma, amp):
    return amp * np.exp(-((grid_x - x0) ** 2 + (grid_y - y0) ** 2) / (2.0 * sigma ** 2))


def simulate_frames(out_dir, seeing_arcsec=1.5, n_stars=3, n_frames=500, config=None,
                    frame_size=(200, 200), prism_sep_pix=40.0, fwhm_pix=3.5,
                    amplitude=15000.0, sky_bg=300.0, read_noise=10.0,
                    altitude_deg=90.0, exptime=0.005, common_motion_pix=3.0,
                    as_cube=False, seed=None, logger=None):
    """Write ``n_frames`` synthetic DIMM frames with a known injected seeing.

    Parameters
    ----------
    out_dir : str
        Directory to write ``f_000000.fits`` ... into (created if needed).
    seeing_arcsec : float
        Injected seeing FWHM (at the observed altitude).
    n_stars : int
        Number of prism doublets placed in the field.
    as_cube : bool
        If True, write a single ``dimm_cube.fits`` instead of per-frame files.
    """
    rng = np.random.default_rng(seed)
    cfg = config or load_config()
    hw = cfg.hardware
    os.makedirs(out_dir, exist_ok=True)

    # Injected seeing -> r0 -> expected differential variance (rad^2 -> pixels).
    seeing_rad = seeing_arcsec * _ARCSEC_TO_RAD
    r0 = _SEEING_CONST * hw.wavelength_m / seeing_rad
    c_l, c_t = response_coefficients(hw.aperture_diameter_m, hw.hole_separation_m)
    var_l_rad2 = 2.0 * hw.wavelength_m ** 2 * r0 ** (-5 / 3) * c_l
    var_t_rad2 = 2.0 * hw.wavelength_m ** 2 * r0 ** (-5 / 3) * c_t
    sig_l_pix = np.sqrt(var_l_rad2) / hw.scale_rad
    sig_t_pix = np.sqrt(var_t_rad2) / hw.scale_rad

    ny, nx = frame_size
    gy, gx = np.mgrid[0:ny, 0:nx]
    sigma_psf = fwhm_pix / 2.355

    # Prism vector horizontal; longitudinal = x, transverse = y.
    pvec = np.array([prism_sep_pix, 0.0])

    # Place well-separated star bases so both spots stay inside the frame.
    margin = 15
    bases = []
    attempts = 0
    while len(bases) < n_stars and attempts < 2000:
        attempts += 1
        bx = rng.uniform(margin, nx - margin - prism_sep_pix)
        by = rng.uniform(margin, ny - margin)
        if all(np.hypot(bx - ox, by - oy) > 1.8 * prism_sep_pix for ox, oy in bases):
            bases.append((bx, by))
    bases = np.array(bases)

    if logger:
        logger.info(f"Simulating {n_frames} frames, {len(bases)} doublets, "
                    f"seeing={seeing_arcsec}\" (r0={r0*100:.1f} cm), "
                    f"sig_l={sig_l_pix:.3f}px sig_t={sig_t_pix:.3f}px")

    t0 = datetime(2026, 8, 27, 20, 0, 0)
    cube = np.zeros((n_frames, ny, nx), dtype=np.float32) if as_cube else None

    for i in range(n_frames):
        clean = np.full((ny, nx), sky_bg, dtype=np.float64)
        for (bx, by) in bases:
            common = rng.normal(0, common_motion_pix, size=2)          # cancels differentially
            dl, dt = rng.normal(0, sig_l_pix), rng.normal(0, sig_t_pix)
            delta = np.array([dl, dt])                                  # long=x, tran=y
            c1 = np.array([bx, by]) + common + 0.5 * delta
            c2 = np.array([bx, by]) + pvec + common - 0.5 * delta
            clean += _render_gaussian(gy, gx, c1[0], c1[1], sigma_psf, amplitude)
            clean += _render_gaussian(gy, gx, c2[0], c2[1], sigma_psf, amplitude)

        noisy = (rng.poisson(np.clip(clean, 0, None))
                 + rng.normal(0, read_noise, (ny, nx))).astype(np.float32)

        if as_cube:
            cube[i] = noisy
        else:
            hdr = _frame_header(t0, i, exptime, altitude_deg, seeing_arcsec, n_stars)
            fits.PrimaryHDU(noisy, hdr).writeto(os.path.join(out_dir, f"f_{i:06d}.fits"),
                                               overwrite=True)

    if as_cube:
        hdr = _frame_header(t0, 0, exptime, altitude_deg, seeing_arcsec, n_stars)
        hdr["FRAMES"] = n_frames
        path = os.path.join(out_dir, "dimm_cube.fits")
        fits.PrimaryHDU(cube, hdr).writeto(path, overwrite=True)
        if logger:
            logger.info(f"Wrote cube {path}")
        return path

    if logger:
        logger.info(f"Wrote {n_frames} frames to {out_dir}")
    return out_dir


def _frame_header(t0, i, exptime, altitude_deg, seeing, n_stars):
    hdr = fits.Header()
    t = t0 + timedelta(seconds=i * (exptime + 0.01))
    hdr["DATE-OBS"] = t.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
    hdr["EXPTIME"] = exptime
    hdr["INSTRUME"] = "SIM_DIMM"
    hdr["OBJCTRA"] = 300.0
    hdr["OBJCTDEC"] = 40.0
    hdr["OBJCTALT"] = altitude_deg
    hdr["SIMSEE"] = (seeing, "Injected seeing [arcsec]")
    hdr["SIMNSTAR"] = (n_stars, "Injected doublet count")
    return hdr
