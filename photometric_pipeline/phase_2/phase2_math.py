# phase2_math.py
import numpy as np
from photutils.background import Background2D, MedianBackground

class MathEngine:
    @staticmethod
    def extract_2d_background(data):
        try:
            bkg = Background2D(data, (50, 50), filter_size=(3, 3), bkg_estimator=MedianBackground())
            return data - bkg.background, bkg.background
        except Exception:
            flat_bg = np.nanmedian(data)
            return data - flat_bg, flat_bg

    @staticmethod
    def calc_scale(target_data, ref_data, target_coords, ref_coords):
        src_x = np.clip(np.round(target_coords[:, 0]).astype(int), 0, target_data.shape[1]-1)
        src_y = np.clip(np.round(target_coords[:, 1]).astype(int), 0, target_data.shape[0]-1)
        dst_x = np.clip(np.round(ref_coords[:, 0]).astype(int), 0, ref_data.shape[1]-1)
        dst_y = np.clip(np.round(ref_coords[:, 1]).astype(int), 0, ref_data.shape[0]-1)
        
        t_flux, r_flux = target_data[src_y, src_x], ref_data[dst_y, dst_x]
        valid = (t_flux > 0) & (r_flux > 0)
        return np.nanmedian(r_flux[valid] / t_flux[valid]) if np.any(valid) else 1.0

    @staticmethod
    def center_crop(data, fraction=0.5):
        """Extracts the center portion of the image array for accurate telemetry/noise analysis."""
        h, w = data.shape
        y1, y2 = int(h * (1 - fraction) / 2), int(h * (1 + fraction) / 2)
        x1, x2 = int(w * (1 - fraction) / 2), int(w * (1 + fraction) / 2)
        return data[y1:y2, x1:x2]