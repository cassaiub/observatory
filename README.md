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

**Linux, macOS and Windows.** One script per platform; it picks an environment,
installs every dependency **including a working plate solver**, installs the
package, registers the Jupyter kernel, and then runs `cassa-doctor` to prove it
worked. Running it again updates in place.

| Platform | Route | Solver you get |
|---|---|---|
| **Linux x86-64** | `./install.sh` | ASTAP, else `solve-field`, else in-process |
| **Linux aarch64** (Raspberry Pi, ARM servers) | `./install.sh --conda` | **ASTAP only** — nothing else is published for ARM Linux |
| **macOS Intel** | `./install.sh` | ASTAP, else `solve-field`, else in-process |
| **macOS Apple Silicon** | `./install.sh` | ASTAP, else in-process (`solve-field` has no `osx-arm64` build) |
| **Windows x64 / ARM64** | `.\install.ps1` | **ASTAP only** — it is the one backend published for Windows |
| **Windows**, alternatively | WSL, then the Linux route | as Linux x86-64 |

**Linux and macOS**

```bash
git clone https://github.com/cassaiub/observatory.git
cd observatory/photometric_pipeline
./install.sh
conda activate cassa-photometry
cassa-doctor
```

**Windows**

```powershell
git clone https://github.com/cassaiub/observatory.git
cd observatory\photometric_pipeline
.\install.ps1
cassa-doctor
```

`install.ps1` does the same things as `install.sh`, with the same options in
PowerShell form (`-Conda`, `-Venv`, `-NewEnv`, `-EnvName`, `-NoDev`,
`-NoDoctor`). WSL still works if you prefer it — it is real x86-64 Linux, so
every instruction applies unchanged inside it, and `.\install.ps1 -Wsl` prints
the setup steps without installing anything.

Install [Miniforge](https://conda-forge.org/download/) first if you have no
conda. **Part II of the documentation** has the step-by-step procedure for each
platform, the macOS Gatekeeper note, and the HPC network-drive workaround.

The other two packages install into the same environment:

```bash
pip install -e DIMM
pip install -e camera_characterization
```

### There is nothing to download for plate solving

Three backends can supply it — **ASTAP**, Astrometry.net's `solve-field`, and an
in-process PyPI solver — and the installer picks the one your platform can run,
trying them in that order. ASTAP leads because it is the only one published for
every platform here: conda-forge builds `astrometry` for `linux-64`/`osx-64`
only, and PyPI's `astrometry` ships no aarch64 and no Windows wheel, so without
ASTAP an ARM Linux box and any Windows machine have **no solver at all**.

Whichever lands, the sky data is fetched **per field** — about 6 MB for ASTAP,
~246 MB for Astrometry.net — and cached under `~/.cache/cassa-photometry/`. An
existing local set is used in preference (`CASSA_ASTROMETRY_INDEX`,
`CASSA_ASTAP_DB`) and nothing is downloaded.

**The products do not depend on which backend solved.** The astrometric residual
`ASTRMS` is measured by the pipeline rather than taken from the solver, and the
header records which route was used as `ASTRMSRC`.

### Choosing how much of the machine to use

Phase 2 holds a whole stack plus its variance planes in memory at once. By
default it asks on a terminal and takes half the cores when there is nobody to
ask. Three settings override that, on the command line or in a config file —
and setting any of them also switches off the prompt:

```bash
cassa-integrate work/phase1 --cpu-fraction 0.5   # half this machine's cores
cassa-integrate work/phase1 --cores 4            # a hard ceiling
cassa-integrate work/phase1 --max-memory-gb 8    # warn before it does not fit
```

### Jupyter

Both installers register a kernel called **`Python (CASSA photometry)`**.
Check that name appears in the top-right of every notebook; if not, *Kernel →
Change Kernel*. Installing packages sets up an environment, it does not make an
existing Jupyter offer it — and the resulting `ModuleNotFoundError: No module
named 'cassa_photometry'` looks exactly like a failed install.

### If something is wrong

```bash
cassa-doctor
```

One line per check — platform, versions, which solver backends are usable and
which one a run would pick, both sky-data caches, network reachability, write
permissions — with the fix for anything that failed. Paste its output when
asking for help.

Large data (raw frames, calibration campaigns, sky databases) is git-ignored;
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
| [`photometric_pipeline/docs/WINDOWS-TESTING.md`](photometric_pipeline/docs/WINDOWS-TESTING.md) | Windows: what is verified, the two platform bugs the first Windows run exposed, and the details worth knowing. |
| [`photometric_pipeline/workshop/`](photometric_pipeline/workshop/) | A lecture, a participant handbook, and five notebooks that reduce a night with known truth. |

## License

MIT — see [`photometric_pipeline/LICENSE`](photometric_pipeline/LICENSE); each
package declares the MIT license in its `pyproject.toml`.
