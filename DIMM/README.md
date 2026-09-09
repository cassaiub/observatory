# cassa-dimm

A professional-grade **Differential Image Motion Monitor (DIMM)** for the CASSA
Observatory. It continuously watches a folder of incoming frames from the
acquisition server, measures the differential motion of each star's two
prism-mask spots, and reports **airmass-corrected atmospheric seeing** as a live
time series.

Hardware target: 8-inch SkyWatcher on an EQ6R-Pro, a two-hole DIMM mask with a
wedge prism over one hole (so every star produces two spots).

## What makes it professional grade

- **Multiple stars.** Detects all spots on a stacked reference, auto-calibrates
  the prism doublet vector from the data, and pairs every star into a doublet —
  then SNR-weight-combines all valid doublets for better statistics (with the
  star-to-star scatter reported as a quality metric).
- **Correct physics.** Sarazin & Roddier (1990) / Tokovinin (2002) reduction
  with the proper aperture term `D^{-1/3}` and separate longitudinal/transverse
  coefficients; longitudinal axis follows the actual prism direction.
- **Airmass correction.** Target altitude from the header (or computed from
  RA/Dec + site) → seeing corrected to the zenith (`eps0 = eps * X^{-0.6}`).
- **Continuous monitoring.** A rolling-window daemon watches the folder, emits an
  estimate every *cadence* frames, and writes a CSV/JSONL log, a live plot, and a
  `status_latest.json` an observatory dashboard can poll.
- **Quality control.** Per-frame sigma-clipping, per-doublet validity fractions,
  longitudinal/transverse consistency, minimum star/frame counts — every window
  is flagged with the reasons it passed or failed.
- **Wide-field ready (0.5°).** Stamp buffering processes each frame on arrival and
  keeps only a single reference plus per-doublet offsets, so memory stays ~one
  frame regardless of window length or field size (streaming 180×1500² frames
  grew RAM ~0.1 GB vs ~1.35 GB for a full-frame buffer). Ellipticity rejection
  drops coma-smeared off-axis spots, and an optional central-field radius confines
  the measurement to the best-corrected part of a wide field.
- **Reproducible & testable.** Central YAML config, a physically-calibrated
  simulator that injects a known seeing, and a test suite that checks the
  reducer recovers it.

## Installation

```bash
conda activate cassa-photometry      # provides numpy/scipy/astropy/photutils/matplotlib
cd DIMM
pip install -e . --no-deps --no-build-isolation
```

This registers three commands: `cassa-dimm-monitor`, `cassa-dimm-batch`,
`cassa-dimm-sim`.

## Usage

```bash
# 1. Copy and edit the config (set your site coordinates + watch folder)
cp config.example.yaml config.yaml

# 2. Run the continuous monitor (the acquisition server writes frames into the folder)
cassa-dimm-monitor --config config.yaml --folder /path/to/watch_folder

# 3. One-shot reduction of an existing cube or a folder of frames
cassa-dimm-batch dimm_simulated_cube.fits

# 4. Generate synthetic data with a known seeing (for testing / demos)
cassa-dimm-sim --out /tmp/dimm_sim --seeing 1.5 --nstars 3 --frames 500
```

The monitor writes into `output.log_dir`:

| File | Contents |
|------|----------|
| `seeing_log.csv` / `.jsonl` | time series: timestamp, seeing (zenith & raw), longitudinal/transverse, r0, airmass, n_frames, n_stars, scatter, SNR, flags |
| `status_latest.json` | the latest measurement (for dashboards) |
| `seeing_timeseries.png` | live seeing-vs-time plot |

## Configuration

All parameters live in `cassa_dimm.config` and are documented in
`config.example.yaml` — optics (D, d, wavelength, pixel/focal length), site
coordinates, watch/window behaviour, detection thresholds, QC limits, and output
paths. Any subset is overridable via the YAML file passed with `--config`.

## Development

```bash
pip install -e ".[dev]"
pytest
```
