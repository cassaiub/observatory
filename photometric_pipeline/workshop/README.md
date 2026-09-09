# cassa-photometry Workshop

Materials for a two-part workshop on the `cassa-photometry` pipeline:

- **Lecture** — `slides/cassa_photometry_lecture.tex` (Beamer): the ideas, operations,
  and stages of the pipeline, with optional math asides.
- **Hands-on** — `notebooks/00…06`: run all four phases on simulated CASSA
  8-inch frames whose **true answers are known**, inspecting the SCI/ERR/DQ
  planes and the catalog at each step, and scoring the result against truth.

Phase 2 plate-solves, and whichever backend your machine installed, the sky data
it needs is **fetched automatically for this field** — about 6 MB of star tiles
for ASTAP, or ~165 MB of index files for Astrometry.net — unless a local set is
already available. Phase 3's zero point queries reference catalogs online the
first time and caches them, after which `--offline` works.

Each phase writes into its own directory under `work/`:
`work/phase1` (calibrated) → `phase2` (masters + WCS) → `phase3` (catalogs) →
`phase4` (diagnostics). Re-running a phase replaces its products in place.

```
workshop/
├── slides/     cassa_photometry_lecture.tex   # Beamer lecture
├── notebooks/  00_setup … 04_phase3_photometry.ipynb
├── raw/        the dataset  (not in git -- generate it, see below)
└── work/       pipeline output       (not in git; created when you run)
```

### The dataset

`raw/` holds whatever tree you put there: the notebooks and the pipeline search
it recursively, so a night straight from the acquisition software —

```
raw/20260903/BIAS/untargeted/   20260903T054529.715412_untargeted_Blue_BIAS_....fits
raw/20260903/DARK/untargeted/   ...
raw/20260903/FLAT/untargeted/   ...
raw/20260903/LIGHT/m22/         20260903T043707.171323_m22_Luminance_LIGHT_....fits
```

— is used as-is, and so is the flat directory the simulator writes. Frames are
sorted by their `IMAGETYP`/`FILTER`/`EXPTIME` headers, never by folder name.
Notebook 00 prints what it found, directory by directory.

**Or generate one — there is nothing to download.**

```bash
cd workshop
cassa-simulate --preset workshop --out .
```

That writes `raw/` (65 frames: 15 science in B/V/R, 30 flats, 10 darks, 10 bias)
plus three truth files. It takes a couple of minutes and is **reproducible from
a seed**, so everyone in the room has bit-identical frames.

The data simulate the CASSA 8-inch f/5 + IMX585 on NGC 7331: one night, a
2048×1400 window (20.4′ × 14.0′) at 0.598″/px. The star field is **real** — Gaia
DR3 positions with synthetic Johnson-Cousins magnitudes — so the frames
genuinely plate-solve and the reference-catalog cross-match genuinely finds the
stars.

What makes it worth learning on is the truth files:

| File | Contents |
|---|---|
| `truth_sources.csv` | every source's position, type and **true magnitude** per band |
| `truth_frames.csv` | every frame's true seeing, transparency, airmass and zero point |
| `truth_config.yaml` | the detector and photometric model used |

So every measurement in the notebooks can be scored, not just inspected. That is
the point: a reduction that *looks* right and a reduction that *is* right are
not distinguishable by eye.

Working with your own data instead? Point `RAW_DIR` at it and skip the
generation step; everything else works the same, minus the truth comparisons.

## Prerequisites (once, per HPC account)

1. **Environment + package** — install [Miniforge](https://conda-forge.org/download/),
   then one command from a clone:
   ```bash
   cd /path/to/photometric_pipeline
   ./install.sh          # Linux, macOS, WSL
   ```
   ```powershell
   cd \path\to\photometric_pipeline
   .\install.ps1         # Windows
   ```
   This brings JupyterLab too, and installs a plate solver: ASTAP first, then
   `solve-field`, then the in-process solver, whichever your platform can run.
   `workshop/handbook/participant_handbook.pdf` walks through it step by step
   for all three platforms; Part II of the observatory manual has the full
   version, including the macOS Gatekeeper note and the HPC workaround.

   The installer also registers the Jupyter kernel. Check the kernel name in the
   top-right of every notebook reads **`Python (CASSA photometry)`**; if not,
   *Kernel → Change Kernel*. A notebook on the wrong kernel reports
   `ModuleNotFoundError: No module named 'cassa_photometry'`, which looks like a
   failed install and is not one.
2. **Star databases / index files** — nothing to do. The pipeline works out what
   this field needs and fetches only that. On a shared machine that already
   holds a full set, point at it instead and nothing is downloaded:
   ```bash
   export CASSA_ASTROMETRY_INDEX=/path/to/astrometry_data   # Astrometry.net
   export CASSA_ASTAP_DB=/path/to/astap_database            # or ASTAP
   ```
3. **Network** — needed once, for notebook `04`'s reference-catalog query and
   for the solver's sky data. Both are cached, so a second pass runs offline. On
   a compute node with no outbound access, run notebooks `03` and `04` once
   somewhere that has it, and use `--offline` afterwards.
4. **Check it all works** before the session:
   ```bash
   cassa-doctor
   ```

## Building the lecture PDF

```bash
cd slides
pdflatex cassa_photometry_lecture.tex   # run twice for the outline/TOC
```

Requires a TeX install with `beamer` (e.g. `texlive-latex-recommended`). Uses only
stock Beamer themes so it compiles anywhere.

## Running the hands-on notebooks

Start Jupyter in the `cassa-photometry` env and run the notebooks **in order**.
Notebook `00` checks the environment and `01` reads only the raw frames; `02`–`04`
each depend on the previous notebook's outputs.

Notebooks `02`-`04` each end with **two exercises that have worked solutions**.
Each asks for one or two specific values, with the answer from the reference
reduction printed beside it, so you can tell at a glance whether your run agrees.

| Notebook | What it does |
|----------|--------------|
| `00_setup` | Verify env, package, plate solver, sky-data cache, dataset. |
| `01_raw_frame_health` | Inventory the raw frames, then health-check them: a go/no-go verdict before anything is reduced. |
| `02_phase1_calibration` | Run ISR → `calibrated_*.fits`; watch ERR/DQ populate, then audit the bad-pixel mask. |
| `03_phase2_integration` | Align + stack + WCS-solve → `Master_*.fits`. Shows the resource budget first: what the run needs, and what share of your machine that is. |
| `04_phase3_photometry` | Zero point, aperture correction, classification → catalog CSV + fluxcal; ends by measuring **your own assigned star**. |

**Before running:** every notebook opens with

```python
from workshop_config import *
```

which resolves `RAW_DIR`, `WORK_DIR` and the per-phase directories relative to
`workshop/`, and fails immediately with a clear message if the dataset has not
been generated. It is one file (`notebooks/workshop_config.py`) rather than the
same thirty lines copied into each notebook, so changing a path changes it
everywhere.

> The dataset is B/V/R. Phases 1–3 handle all three filters in one pass; the
> notebooks inspect one filter at a time and say which.
