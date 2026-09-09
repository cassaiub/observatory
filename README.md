# CASSA Observatory

Software for the **CASSA Observatory** — an 8-inch SkyWatcher on an EQ6R-Pro,
with an iTelescope-compatible imaging workflow, a DIMM seeing monitor, and a
CMOS sensor-characterization bench. This repository collects three independent,
installable Python packages plus a shared master documentation set.

## Packages

| Package | Folder | Commands | What it does |
|---------|--------|----------|--------------|
| **cassa-photometry** | [`photometric_pipeline/`](photometric_pipeline/README.md) | `cassa-calibrate`, `cassa-integrate`, `cassa-photometry`, `cassa-diagnose`, `cassa-verify`, `cassa-run` | End-to-end imaging reduction (Phases 1–4): calibration → stacking + WCS → photometry/zero point → diagnostics, with a full `SCI`/`ERR`/`DQ` error budget and a configurable step plan (exclude, reorder or extend any phase's steps). |
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

**Supported platforms: Linux and macOS**, on both Intel and ARM. Windows works
through WSL, which is real x86-64 Linux — every instruction then applies
unchanged inside it.

The photometric pipeline installs itself:

```bash
git clone https://github.com/cassaiub/observatory.git
cd observatory/photometric_pipeline
./install.sh
conda activate cassa-photometry
cassa-doctor
```

`install.sh` picks an environment, installs every dependency **including a
working plate solver**, installs the package, and then runs `cassa-doctor` to
prove it worked. Part II of the documentation has the step-by-step procedure for
each of the three operating systems, plus the HPC network-drive workaround.

The other two packages install into the same environment:

```bash
pip install -e DIMM
pip install -e camera_characterization
```

Notes:
- **There is nothing to download for plate solving.** Three backends can supply
  it — ASTAP, Astrometry.net's `solve-field`, and an in-process PyPI solver —
  and the installer picks the one your platform can run. Whichever it is, the
  sky data is fetched **per field**: about 6 MB for ASTAP, ~246 MB for
  Astrometry.net, cached under `~/.cache/cassa-photometry/`. An existing local
  set is used in preference (`CASSA_ASTROMETRY_INDEX`, `CASSA_ASTAP_DB`).
- The reduction **products do not depend on which backend solved**: the
  astrometric residual `ASTRMS` is measured by the pipeline rather than taken
  from the solver, and the header records which route was used as `ASTRMSRC`.
- Large data (raw frames, calibration campaigns, sky databases) is git-ignored;
  only source, configs, and documentation are tracked.

## Documentation

The full **CASSA Observatory Master Documentation** lives in [`docs/`](docs/)
(`main.tex`, compiled to `main.pdf`) and is organised as:

- **Part I** — Astronomical CMOS Sensor Characterization SOP (capture procedures)
- **Part II** — Pipeline Installation & Deployment, for Linux, macOS and Windows
- **Part III** — Image Processing Pipeline Architecture (Phases 1–4)
- **Part IV** — Atmospheric Seeing Monitor (`cassa-dimm`)
- **Part V** — Sensor Characterization Pipeline (`cassa-camchar`)
- **Part VI** — Complete Reference: every command and flag, every configuration
  key with its default, every environment variable, every FITS keyword and
  catalog column, and the module-by-module Python API

Rebuild it with `latexmk -pdf main.tex`.

Shorter, task-focused documents live with the pipeline itself:

| Document | What it covers |
|---|---|
| [`photometric_pipeline/README.md`](photometric_pipeline/README.md) | What the pipeline is, and how to run it. |
| [`photometric_pipeline/docs/REFERENCE.md`](photometric_pipeline/docs/REFERENCE.md) | The same reference as Part VI, in markdown. |
| [`photometric_pipeline/docs/CUSTOMIZING.md`](photometric_pipeline/docs/CUSTOMIZING.md) | Excluding, reordering and adding reduction steps; describing your own setup; changing an algorithm and staying mergeable with upstream. |
| [`photometric_pipeline/docs/CHANGELOG.md`](photometric_pipeline/docs/CHANGELOG.md) | What changed, and what it means for existing data. |
| [`photometric_pipeline/workshop/`](photometric_pipeline/workshop/) | A lecture, a participant handbook, and five notebooks that reduce a night with known truth. |

## License

MIT — see [`photometric_pipeline/LICENSE`](photometric_pipeline/LICENSE); each
package declares the MIT license in its `pyproject.toml`.
