"""Watch-folder seeing monitor (continuous) and one-shot batch reduction.

The monitor polls a folder for newly-arrived individual FITS frames, keeps a
rolling window of the most recent frames, and every ``cadence`` frames emits a
seeing estimate to a CSV/JSONL log, a live plot, and a ``status_latest.json``
that an observatory dashboard can poll. Uses stdlib polling (no watchdog dep).
"""

import os
import csv
import json
import glob
import time
from dataclasses import asdict
from datetime import datetime, timezone

from cassa_dimm.config import load_config
from cassa_dimm.logging_utils import get_logger
from cassa_dimm.io import read_frame, read_cube, extract_meta, file_is_stable
from cassa_dimm.stream import StreamEstimator
from cassa_dimm.window import estimate_window

_CSV_FIELDS = ["timestamp", "seeing_zenith", "seeing_raw", "seeing_l", "seeing_t",
               "r0_cm", "airmass", "n_frames", "n_stars", "star_scatter",
               "mean_snr", "seeing_err", "flags"]


class DimmMonitor:
    """Continuous rolling-window seeing monitor over a watched folder."""

    def __init__(self, config=None, logger=None):
        self.config = config or load_config()
        self.log = logger or get_logger("cassa_dimm", log_dir=self.config.output.log_dir,
                                        level=self.config.log_level)
        # Stamp-buffering estimator: memory ~= one frame + a few floats per doublet,
        # independent of window length or field size (no full-frame window buffer).
        self.estimator = StreamEstimator(self.config, logger=self.log)
        self.seen = set()
        self.new_since_emit = 0
        self.hist_times, self.hist_seeing, self.hist_err = [], [], []

        out = self.config.output
        os.makedirs(out.log_dir, exist_ok=True)
        self.csv_path = os.path.join(out.log_dir, out.csv_name)
        self.jsonl_path = os.path.join(out.log_dir, out.jsonl_name)
        self.status_path = os.path.join(out.log_dir, out.status_name)
        self.plot_path = os.path.join(out.log_dir, out.plot_name)
        if not os.path.exists(self.csv_path):
            with open(self.csv_path, "w", newline="") as fh:
                csv.DictWriter(fh, fieldnames=_CSV_FIELDS).writeheader()

    # -- main loop ------------------------------------------------------------
    def run(self):
        w = self.config.watch
        folder = os.path.abspath(w.folder)
        os.makedirs(folder, exist_ok=True)
        self.log.info(f"DIMM monitor watching {folder} (window={w.window_frames}, "
                      f"cadence={w.cadence_frames})")

        if not w.process_existing:
            self.seen.update(glob.glob(os.path.join(folder, w.pattern)))
            self.log.info(f"Ignoring {len(self.seen)} pre-existing files.")

        try:
            while True:
                self._scan(folder)
                time.sleep(w.poll_interval_s)
        except KeyboardInterrupt:
            self.log.info("Monitor stopped by user.")

    def _scan(self, folder):
        w = self.config.watch
        new_files = sorted(f for f in glob.glob(os.path.join(folder, w.pattern))
                          if f not in self.seen)
        for path in new_files:
            if not file_is_stable(path, w.stabilization_s):
                continue  # still being written; pick it up next scan
            self.seen.add(path)
            if self._ingest(path):
                self.new_since_emit += 1

        if (self.new_since_emit >= w.cadence_frames
                and self.estimator.max_track_len() >= self.config.qc.min_frames):
            self._emit()
            self.new_since_emit = 0

    def _ingest(self, path):
        try:
            data, header = read_frame(path)
        except Exception as exc:
            self.log.warning(f"Skipping {os.path.basename(path)}: {exc}")
            return False
        self.estimator.add_frame(data, extract_meta(header))
        return True

    # -- emit an estimate -----------------------------------------------------
    def _emit(self):
        result = self.estimator.estimate()
        ts = self.estimator.latest_time() or datetime.now(timezone.utc).isoformat()
        record = self._record(ts, result)

        self._append_csv(record)
        self._append_jsonl(record)
        self._write_status(record)
        self._update_plot(ts, result)

        if result.n_stars > 0:
            self.log.info(f"[{ts}] seeing={result.seeing_zenith:.2f}\" "
                          f"+/-{result.seeing_err:.2f} | stars={result.n_stars} "
                          f"| airmass={result.airmass:.2f} | flags={result.flags}")
        else:
            self.log.warning(f"[{ts}] no valid measurement | flags={result.flags}")

    @staticmethod
    def _record(ts, result):
        rec = asdict(result)
        rec["timestamp"] = ts
        rec["flags"] = ";".join(result.flags)
        return {k: rec.get(k) for k in _CSV_FIELDS}

    def _append_csv(self, record):
        with open(self.csv_path, "a", newline="") as fh:
            csv.DictWriter(fh, fieldnames=_CSV_FIELDS).writerow(record)

    def _append_jsonl(self, record):
        with open(self.jsonl_path, "a") as fh:
            fh.write(json.dumps(_json_safe(record)) + "\n")

    def _write_status(self, record):
        status = _json_safe(record)
        status["updated"] = datetime.now(timezone.utc).isoformat()
        with open(self.status_path, "w") as fh:
            json.dump(status, fh, indent=2)

    def _update_plot(self, ts, result):
        if result.n_stars == 0:
            return
        from cassa_dimm.plots import plot_timeseries
        try:
            self.hist_times.append(datetime.fromisoformat(ts.replace("Z", "+00:00")))
        except (ValueError, AttributeError):
            self.hist_times.append(datetime.now(timezone.utc))
        self.hist_seeing.append(result.seeing_zenith)
        self.hist_err.append(result.seeing_err)
        plot_timeseries(self.hist_times, self.hist_seeing, self.hist_err,
                        self.plot_path, max_points=self.config.output.plot_max_points)


def _json_safe(record):
    out = {}
    for k, v in record.items():
        if isinstance(v, float) and (v != v):  # NaN
            out[k] = None
        else:
            out[k] = v
    return out


# -- one-shot batch reduction ------------------------------------------------- #
def run_batch(path, config=None, logger=None):
    """Reduce a single cube file or a directory of frames to one seeing value."""
    config = config or load_config()
    logger = logger or get_logger("cassa_dimm")
    path = os.path.abspath(path)

    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, config.watch.pattern)))
        frames, metas = [], []
        for f in files:
            try:
                data, header = read_frame(f)
            except Exception:
                continue
            frames.append(data)
            metas.append(extract_meta(header))
    else:
        cube, header = read_cube(path)
        frames = list(cube)
        metas = [extract_meta(header)] * len(frames)

    logger.info(f"Batch: {len(frames)} frames from {path}")
    result = estimate_window(frames, metas, config, logger)
    logger.info(f"Seeing (zenith): {result.seeing_zenith:.2f}\" +/- {result.seeing_err:.2f} "
                f"| raw={result.seeing_raw:.2f}\" | r0={result.r0_cm:.1f} cm "
                f"| stars={result.n_stars} | flags={result.flags}")
    return result
