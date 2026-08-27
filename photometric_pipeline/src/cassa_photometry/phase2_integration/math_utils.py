"""Numerical helpers for the integration phase.

Adds error-aware routines on top of the original background/scale utilities:
plane registration (used for both SCI and variance) and an inverse-variance
weighted stack that returns a propagated variance plane.
"""

import numpy as np
from astropy.stats import sigma_clip
from skimage.transform import warp
from photutils.background import Background2D, MedianBackground


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
    def calc_scale(target_data, ref_data, target_coords, ref_coords):
        """Median flux ratio (reference / target) at matched star positions."""
        src_x = np.clip(np.round(target_coords[:, 0]).astype(int), 0, target_data.shape[1] - 1)
        src_y = np.clip(np.round(target_coords[:, 1]).astype(int), 0, target_data.shape[0] - 1)
        dst_x = np.clip(np.round(ref_coords[:, 0]).astype(int), 0, ref_data.shape[1] - 1)
        dst_y = np.clip(np.round(ref_coords[:, 1]).astype(int), 0, ref_data.shape[0] - 1)

        t_flux, r_flux = target_data[src_y, src_x], ref_data[dst_y, dst_x]
        valid = (t_flux > 0) & (r_flux > 0)
        return np.nanmedian(r_flux[valid] / t_flux[valid]) if np.any(valid) else 1.0

    @staticmethod
    def register_plane(plane, inverse_map, output_shape, order=3, cval=np.nan):
        """Warp a 2D plane onto the reference grid."""
        return warp(plane, inverse_map=inverse_map, output_shape=output_shape,
                    order=order, cval=cval, preserve_range=True)

    @staticmethod
    def register_variance(variance, inverse_map, output_shape, cval=np.nan):
        """Warp a variance plane (bilinear, clipped non-negative to avoid overshoot)."""
        warped = warp(variance, inverse_map=inverse_map, output_shape=output_shape,
                      order=1, cval=cval, preserve_range=True)
        return np.clip(warped, 0.0, None)

    @staticmethod
    def weighted_stack(cube, varcube, weights, sigma=3.0, maxiters=3, use_sigma_clip=False):
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
            If True, reject outliers with sigma-clipping; otherwise reject the
            per-pixel min and max (classic min/max rejection for small stacks).

        Returns
        -------
        master : ndarray, shape (H, W)
            Weighted mean; NaN where every frame was rejected.
        master_var : ndarray, shape (H, W)
            Propagated variance of the weighted mean:
            ``sum(w_i^2 var_i) / (sum w_i)^2`` over the surviving frames.
        """
        cube = np.asarray(cube, dtype=np.float64)
        varcube = np.asarray(varcube, dtype=np.float64)
        w = np.asarray(weights, dtype=np.float64).reshape(-1, 1, 1)

        if use_sigma_clip:
            clipped = sigma_clip(cube, sigma=sigma, maxiters=maxiters, axis=0, masked=True)
            reject = np.ma.getmaskarray(clipped)
        else:
            reject = (cube == np.nanmin(cube, axis=0)) | (cube == np.nanmax(cube, axis=0))
        valid = reject == False  # noqa: E712 -- explicit boolean array
        valid &= np.isfinite(cube) & np.isfinite(varcube)

        w_valid = np.where(valid, w, 0.0)
        wsum = w_valid.sum(axis=0)
        data_num = np.where(valid, w * cube, 0.0).sum(axis=0)
        var_num = np.where(valid, (w ** 2) * varcube, 0.0).sum(axis=0)

        with np.errstate(divide="ignore", invalid="ignore"):
            master = np.where(wsum > 0, data_num / wsum, np.nan)
            master_var = np.where(wsum > 0, var_num / (wsum ** 2), np.nan)
        return master, master_var

    @staticmethod
    def center_crop(data, fraction=0.5):
        """Extract the central portion of an array for telemetry/noise analysis."""
        h, w = data.shape
        y1, y2 = int(h * (1 - fraction) / 2), int(h * (1 + fraction) / 2)
        x1, x2 = int(w * (1 - fraction) / 2), int(w * (1 + fraction) / 2)
        return data[y1:y2, x1:x2]
