# CASSA Observatory

Software for the **CASSA Observatory** — an 8-inch SkyWatcher on an EQ6R-Pro,
with an iTelescope-compatible imaging workflow, a DIMM seeing monitor, and a
CMOS sensor-characterization bench. This repository collects three independent,
installable Python packages plus a shared master documentation set.

## Packages

| Package | Folder | Commands | What it does |
|---------|--------|----------|--------------|
| **cassa-photometry** | [`photometric_pipeline/`](photometric_pipeline/README.md) | `cassa-calibrate`, `cassa-integrate`, `cassa-photometry`, `cassa-diagnose`, `cassa-verify`, `cassa-run` | End-to-end imaging reduction (Phases 1–4): calibration → stacking + WCS → photometry/zero point → diagnostics, with a full `SCI`/`ERR`/`DQ` error budget. |
| **cassa-dimm** | [`DIMM/`](DIMM/README.md) | `cassa-dimm-monitor`, `cassa-dimm-batch`, `cassa-dimm-sim` | Continuous atmospheric-**seeing** monitor: watches a folder of frames, measures the differential motion of prism-mask star doublets, and reports airmass-corrected seeing as a live time series. |
| **cassa-camchar** | [`camera_characterization/`](camera_characterization/README.md) | `cassa-camchar-analyze`, `cassa-camchar-sim` | CMOS **sensor characterization** via the Photon Transfer Curve: gain, read noise, full well, dynamic range, dark current, QE, and filter transmission. |

Each package has its own README with detailed usage; the [`docs/`](docs/) folder
holds the observatory **master documentation**.

## Repository layout

```
cassa_observatory/
├── photometric_pipeline/     # cassa-photometry  — imaging reduction (Phases 1–4)
├── DIMM/                      # cassa-dimm        — atmospheric seeing monitor
├── camera_characterization/  # cassa-camchar     — CMOS sensor characterization
└── docs/                      # master documentation (LaTeX source + compiled PDF)
```

## Installation

All three packages share one Conda environment (scientific libraries + the
`solve-field` binary that the imaging and DIMM pipelines need). Build it once
(see Part II of the documentation for the `environment.yml` and any HPC notes),
then install whichever packages you need in editable mode:

```bash
conda activate pipeline_env

pip install -e photometric_pipeline    --no-deps --no-build-isolation
pip install -e DIMM                    --no-deps --no-build-isolation
pip install -e camera_characterization --no-deps --no-build-isolation
```

Notes:
- **Astrometry.net** (`solve-field`) is a system binary, not a `pip` package —
  install it via conda-forge (`astrometry`) or apt (`astrometry.net`). The
  plate-solving in Phase 2 and any DIMM astrometry rely on it.
- The astrometry **index files** (several GB) are downloaded separately and
  pointed to via `CASSA_ASTROMETRY_INDEX` (see `photometric_pipeline/README.md`).
- Large data (raw frames, calibration campaigns, astrometry indices) is
  git-ignored; only source, configs, and documentation are tracked.

## Documentation

The full **CASSA Observatory Master Documentation** lives in [`docs/`](docs/)
(`main.tex`, compiled to `main.pdf`) and is organised as:

- **Part I** — Astronomical CMOS Sensor Characterization SOP (capture procedures)
- **Part II** — Pipeline Installation & Deployment
- **Part III** — Image Processing Pipeline Architecture (Phases 1–4)
- **Part IV** — Atmospheric Seeing Monitor (`cassa-dimm`)
- **Part V** — Sensor Characterization Pipeline (`cassa-camchar`)

Rebuild it with `pdflatex main.tex` (run twice for cross-references).

## License

MIT — see [`photometric_pipeline/LICENSE`](photometric_pipeline/LICENSE); each
package declares the MIT license in its `pyproject.toml`.
