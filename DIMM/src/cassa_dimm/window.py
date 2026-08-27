"""Batch (full-window) seeing estimate.

Thin wrapper over :class:`~cassa_dimm.stream.StreamEstimator`: it fixes the
doublet geometry once from the full-window mean, then feeds every frame through
the same per-frame centroiding path used by the live monitor.
"""

import numpy as np

from cassa_dimm.seeing import WindowResult
from cassa_dimm.stream import StreamEstimator


def _empty_result(n_frames, flag):
    return WindowResult(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                        n_frames, 0, np.nan, np.nan, np.nan, [flag])


def estimate_window(frames, metas, config, logger=None):
    """Estimate seeing from a list of 2D frames (and optional per-frame metas)."""
    n_frames = len(frames)
    if n_frames < config.qc.min_frames:
        return _empty_result(n_frames, "TOO_FEW_FRAMES")

    est = StreamEstimator(config, logger, max_window=n_frames)
    reference = np.mean(np.asarray(frames, dtype=np.float64), axis=0)
    if not est.set_geometry(reference):
        return _empty_result(n_frames, est.status)
    for d, m in zip(frames, metas):
        est.add_frame(d, m, allow_refresh=False)
    return est.estimate()
