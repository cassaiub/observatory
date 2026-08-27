"""Incremental, memory-bounded DIMM estimator (stamp buffering).

Instead of buffering a window of full frames (which for a 0.5-degree field is
gigabytes of RAM), the estimator processes each frame as it arrives: it keeps a
single exponential-moving-average reference for detection and, per doublet, only
a rolling deque of the differential centroid offsets ``(dx, dy)``. Memory is
therefore ~one frame plus a few floats per doublet, independent of the window
length or field size.

The same class powers batch mode (``estimate_window``), where the geometry is
fixed once from the full-window mean.
"""

from collections import deque

import numpy as np

from cassa_dimm.detect import detect_sources, estimate_prism_vector, pair_doublets
from cassa_dimm.centroid import centroid_in_box
from cassa_dimm.seeing import (
    measure_doublet, combine_window, kasten_young_airmass, WindowResult,
)
from cassa_dimm.io import altitude_deg


class StreamEstimator:
    """Rolling, per-frame DIMM estimator with bounded memory."""

    def __init__(self, config, logger=None, max_window=None):
        self.config = config
        self.logger = logger
        self.max_window = max_window if max_window is not None else config.watch.window_frames
        self.ref = None                 # rolling reference (single array)
        self._ref_n = 0                 # frames averaged so far (for the warm-up mean)
        self.geometry = []              # list of Doublet
        self.prism_unit = None
        self.tracks = []                # list of deque[(dx, dy)], aligned with geometry
        self.metas = deque(maxlen=self.max_window)
        self.status = "OK"
        self.frame_count = 0
        self.frames_since_refresh = 0
        self.lost = 0

    # -- geometry ------------------------------------------------------------
    def set_geometry(self, reference):
        """(Re)detect doublets and the prism vector from a reference image.

        Existing per-doublet tracks are carried over to spatially-matched new
        doublets so the rolling window survives a re-detection.
        """
        sources = detect_sources(reference, self.config)
        if len(sources) < 2:
            self.status = "NO_SOURCES"
            self.geometry, self.tracks, self.prism_unit = [], [], None
            return False
        pv = estimate_prism_vector(sources, self.config)
        if pv is None or float(np.hypot(*pv[0])) <= 0:
            self.status = "NO_PRISM_VECTOR"
            self.geometry, self.tracks, self.prism_unit = [], [], None
            return False
        vec = pv[0]
        norm = float(np.hypot(*vec))
        self.prism_unit = (vec[0] / norm, vec[1] / norm)

        new_doublets = pair_doublets(sources, vec, self.config)
        if not new_doublets:
            self.status = "NO_DOUBLETS"
            self.geometry, self.tracks = [], []
            return False

        old_geo, old_tracks = self.geometry, self.tracks
        new_tracks = []
        for db in new_doublets:
            carried = deque(maxlen=self.max_window)
            for og, ot in zip(old_geo, old_tracks):   # carry history from a near match
                if np.hypot(db.x1 - og.x1, db.y1 - og.y1) < self.config.detection.centroid_box:
                    carried = deque(ot, maxlen=self.max_window)
                    break
            new_tracks.append(carried)

        self.geometry, self.tracks = new_doublets, new_tracks
        self.frames_since_refresh, self.lost, self.status = 0, 0, "OK"
        if self.logger:
            self.logger.debug(f"geometry: {len(self.geometry)} doublets | prism={vec.round(2)}")
        return True

    # -- ingest --------------------------------------------------------------
    def add_frame(self, data, meta, allow_refresh=True):
        """Process one frame: update the reference and append per-doublet offsets."""
        data = np.asarray(data, dtype=np.float32)
        self._update_reference(data)
        self.frame_count += 1
        self.frames_since_refresh += 1

        if allow_refresh and not self.geometry and self.frame_count >= self.config.detection.warmup_frames:
            self.set_geometry(self.ref)

        if self.geometry:
            box = self.config.detection.centroid_box
            success = 0
            for i, db in enumerate(self.geometry):
                c1 = centroid_in_box(data, db.x1, db.y1, box)
                c2 = centroid_in_box(data, db.x2, db.y2, box)
                if c1 is not None and c2 is not None:
                    self.tracks[i].append((c1[0] - c2[0], c1[1] - c2[1]))
                    success += 1
            self.metas.append(meta)
            self.lost = self.lost + 1 if success < 0.5 * len(self.geometry) else 0

            if allow_refresh and (self.lost >= 8
                                  or self.frames_since_refresh >= self.config.detection.refresh_frames):
                self.set_geometry(self.ref)

    def _update_reference(self, data):
        """Running mean for the first ``ema_window`` frames, then an EMA.

        The warm-up mean gives a clean, high-SNR reference for the first
        detection; the subsequent EMA lets it track slow field drift. Both are
        O(1) in memory (a single array).
        """
        window = max(self.config.detection.ema_window, 1)
        if self.ref is None:
            self.ref = data.astype(np.float64)
            self._ref_n = 1
        elif self._ref_n < window:
            self._ref_n += 1
            self.ref += (data - self.ref) / self._ref_n     # incremental mean
        else:
            alpha = 1.0 / window
            self.ref *= (1.0 - alpha)
            self.ref += alpha * data

    # -- estimate ------------------------------------------------------------
    def max_track_len(self):
        return max((len(t) for t in self.tracks), default=0)

    def latest_time(self):
        for m in reversed(self.metas):
            if m is not None and m.time is not None:
                return m.time.isot
        return None

    def estimate(self):
        """Combine all doublets with enough history into a window result."""
        if not self.geometry:
            return WindowResult(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                                self.frame_count, 0, np.nan, np.nan, np.nan, [self.status])
        results = []
        lengths = []
        for i, db in enumerate(self.geometry):
            t = self.tracks[i]
            lengths.append(len(t))
            if len(t) < self.config.qc.min_frames:
                continue
            arr = np.array(t)
            results.append(measure_doublet(arr[:, 0], arr[:, 1], self.prism_unit,
                                           self.config.hardware, self.config.qc, snr=db.snr))
        n_frames = int(np.median(lengths)) if lengths else 0
        return combine_window(results, self._airmass(), n_frames, self.config.qc)

    def _airmass(self):
        if not self.metas:
            return None
        mid = list(self.metas)[len(self.metas) // 2]
        if mid is None:
            return None
        alt = altitude_deg(mid, self.config.site)
        return kasten_young_airmass(alt) if alt is not None else None
