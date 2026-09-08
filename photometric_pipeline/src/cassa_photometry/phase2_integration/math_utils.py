"""Numerical helpers for the integration phase.

Adds error-aware routines on top of the original background/scale utilities:
plane registration (used for both SCI and variance) and an inverse-variance
weighted stack that returns a propagated variance plane.
"""

import numpy as np
from astropy.stats import sigma_clip
from photutils.background import Background2D, MedianBackground
from skimage.transform import warp


class MathEngine:
    @staticmethod
    def extract_2d_background(data, box_size=(50, 50), filter_size=(3, 3)):
        """Return ``(data - background, background)`` using a 2D estimator."""
        try:
            bkg = Background2D(data, box_size, filter_size=filter_size,
                              bkg_estimator=MedianBackground())
            return data - bkg.background, bkg.background
        except Exception:
            flat_bg = np.nanmedian(data)
            return data - flat_bg, flat_bg

    @staticmethod
    def calc_scale(target_data, ref_data, target_coords, ref_coords, aperture_radius=5.0):
        """Transparency ratio (reference / target) from matched stars.

        Measured as a **sum over an aperture**, not from the peak pixel. A star's
        peak scales roughly as ``1/FWHM**2``, so a peak-pixel ratio reads a
        change in seeing as a change in transparency -- and that spurious factor
        then multiplies the whole frame before stacking, carrying straight into
        the master and hence the zero point. Integrated flux is conserved under
        seeing changes, which is the entire point of using it.
        """
        radius = max(float(aperture_radius), 1.5)
        target_flux = _aperture_sums(target_data, target_coords, radius)
        ref_flux = _aperture_sums(ref_data, ref_coords, radius)

        valid = (
            np.isfinite(target_flux) & np.isfinite(ref_flux)
            & (target_flux > 0) & (ref_flux > 0)
        )
        if not np.any(valid):
            return 1.0
        ratios = ref_flux[valid] / target_flux[valid]
        # Sigma-clip: a mismatched pair or a cosmic ray on one star should not
        # rescale the frame.
        if ratios.size >= 4:
            clipped = sigma_clip(ratios, sigma=3.0, maxiters=3, masked=True)
            kept = ratios[~np.ma.getmaskarray(clipped)]
            if kept.size:
                ratios = kept
        return float(np.median(ratios))

    @staticmethod
    def register_plane(plane, inverse_map, output_shape, order=3, cval=np.nan):
        """Warp a 2D plane onto the reference grid."""
        return warp(plane, inverse_map=inverse_map, output_shape=output_shape,
                    order=order, cval=cval, preserve_range=True)

    @staticmethod
    def register_variance(variance, inverse_map, output_shape, order=1, cval=np.nan):
        """Warp a variance plane onto the reference grid.

        Two caveats, stated rather than hidden:

        * Interpolation forms each output pixel as ``sum(w_i x_i)``, whose
          variance is ``sum(w_i**2 var_i)`` -- not ``sum(w_i var_i)``, which is
          what warping the variance plane directly computes. The result is
          therefore an over-estimate of the per-pixel variance.
        * It also **correlates neighbouring pixels**, so summing this plane over
          an aperture under-estimates the true uncertainty of that sum. Phase 2
          records the interpolation order used so a consumer can know.

        Bilinear (order 1) is kept deliberately even though the data is warped
        bicubic: a bicubic kernel rings, and negative variance is meaningless.
        """
        warped = warp(variance, inverse_map=inverse_map, output_shape=output_shape,
                      order=order, cval=cval, preserve_range=True)
        return np.clip(warped, 0.0, None)

    @staticmethod
    def register_mask(mask, inverse_map, output_shape):
        """Warp a boolean/integer mask with nearest-neighbour interpolation.

        Nearest-neighbour because a DQ bit is a statement about a pixel, not a
        quantity to be averaged: interpolating it would invent fractional flags
        and spread a single bad pixel over its neighbours' values.
        """
        warped = warp(np.asarray(mask, dtype=float), inverse_map=inverse_map,
                      output_shape=output_shape, order=0, cval=0.0,
                      preserve_range=True)
        return np.rint(warped).astype(np.int32)

    @staticmethod
    def weighted_stack(cube, varcube, weights, sigma=3.0, maxiters=3,
                       use_sigma_clip=False, bad_pixels=None):
        """Inverse-variance weighted combine with rejection, returning SCI and VAR.

        Parameters
        ----------
        cube : ndarray, shape (N, H, W)
            Aligned, flux-scaled science frames.
        varcube : ndarray, shape (N, H, W)
            Per-pixel variances matching ``cube``.
        weights : array_like, shape (N,)
            Per-frame scalar weights (typically ``1/noise**2``).
        sigma, maxiters : float, int
            Sigma-clip parameters (used when ``use_sigma_clip``).
        use_sigma_clip : bool
            If True, reject outliers with sigma-clipping. If False, nothing is
            rejected: with only a handful of frames there is not enough
            information to identify an outlier without throwing away most of the
            signal.
        bad_pixels : ndarray of bool, shape (N, H, W), optional
            Per-frame mask of pixels to exclude (from the DQ planes).

        Returns
        -------
        master : ndarray, shape (H, W)
            Weighted mean; NaN where every frame was rejected.
        master_var : ndarray, shape (H, W)
            Propagated variance of the weighted mean:
            ``sum(w_i^2 var_i) / (sum w_i)^2`` over the surviving frames.
        n_used : ndarray, shape (H, W)
            How many frames actually contributed to each pixel. Zero means no
            coverage, and is what the master's ``DQ_NO_DATA`` bit is set from.
        """
        cube = np.asarray(cube, dtype=np.float64)
        varcube = np.asarray(varcube, dtype=np.float64)
        w = np.asarray(weights, dtype=np.float64).reshape(-1, 1, 1)

        if use_sigma_clip:
            clipped = sigma_clip(cube, sigma=sigma, maxiters=maxiters, axis=0, masked=True)
            reject = np.ma.getmaskarray(clipped)
        else:
            # No rejection below the sigma-clip threshold. Min/max rejection --
            # what this used to do -- discards two frames of every three, so a
            # 3-frame stack kept ONE frame per pixel, and where several frames
            # shared a value (flat regions after float32 quantisation, saturated
            # cores, zeroed no-data) it could reject them all and leave a NaN.
            reject = np.zeros(cube.shape, dtype=bool)

        valid = ~reject
        valid &= np.isfinite(cube) & np.isfinite(varcube)
        if bad_pixels is not None:
            # Pixels phase 1 flagged are excluded from the combine rather than
            # averaged in as good data.
            valid &= ~np.asarray(bad_pixels, dtype=bool)

        w_valid = np.where(valid, w, 0.0)
        wsum = w_valid.sum(axis=0)
        data_num = np.where(valid, w * cube, 0.0).sum(axis=0)
        var_num = np.where(valid, (w ** 2) * varcube, 0.0).sum(axis=0)

        with np.errstate(divide="ignore", invalid="ignore"):
            master = np.where(wsum > 0, data_num / wsum, np.nan)
            master_var = np.where(wsum > 0, var_num / (wsum ** 2), np.nan)
        return master, master_var, valid.sum(axis=0).astype(np.int32)

    @staticmethod
    def center_crop(data, fraction=0.5):
        """Extract the central portion of an array for telemetry/noise analysis."""
        h, w = data.shape
        y1, y2 = int(h * (1 - fraction) / 2), int(h * (1 + fraction) / 2)
        x1, x2 = int(w * (1 - fraction) / 2), int(w * (1 + fraction) / 2)
        return data[y1:y2, x1:x2]


def _aperture_sums(data, coords, radius):
    """Background-subtracted flux in a circular aperture at each position.

    A small local median annulus is subtracted so the ratio measures starlight
    rather than sky, which differs between frames.
    """
    from photutils.aperture import (
        ApertureStats,
        CircularAnnulus,
        CircularAperture,
        aperture_photometry,
    )

    positions = np.asarray(coords, dtype=float)
    if positions.size == 0:
        return np.array([])

    aperture = CircularAperture(positions, r=radius)
    annulus = CircularAnnulus(positions, r_in=radius * 2.0, r_out=radius * 3.0)
    clean = np.nan_to_num(np.asarray(data, dtype=float), nan=0.0)
    try:
        table = aperture_photometry(clean, aperture)
        background = ApertureStats(clean, annulus).median
        return np.asarray(table["aperture_sum"]) - background * aperture.area
    except Exception:
        return np.full(len(positions), np.nan)
