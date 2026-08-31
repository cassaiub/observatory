# cassa-photometry Workshop

Materials for a two-part workshop on the `cassa-photometry` pipeline:

- **Lecture** — `slides/cassa_photometry_lecture.tex` (Beamer): the ideas, operations,
  and stages of the pipeline, with optional math asides.
- **Hands-on** — `notebooks/00…05`: run all four phases on real frames, inspecting
  the SCI/ERR/DQ planes and the catalog at each step.

Phase 2 plate-solves against a **shared local Astrometry.net index** (no network).
Phase 3's zero-point cross-match queries online reference catalogs
(APASS/Pan-STARRS/SDSS), so that one step needs outbound internet on the node.

Each phase writes into its own directory under `work/`:
`work/phase1` (calibrated) → `phase2` (masters + WCS) → `phase3` (catalogs) →
`phase4` (diagnostics). Re-running a phase replaces its products in place.

```
workshop/
├── slides/     cassa_photometry_lecture.tex   # Beamer lecture
├── notebooks/  00_setup … 05_phase4_diagnostics.ipynb
├── raw/        the provided dataset  (not in git -- see below)
└── work/       pipeline output       (not in git; created when you run)
```

### The dataset

`raw/` is **not committed** (155 MB of FITS). Obtain the workshop frames
separately and unpack them into `workshop/raw/`, or point `RAW_DIR` at wherever
you keep them. The set is one night of NGC 7331 in B/V/R plus its calibration
frames: 14 science, 43 flats, 10 darks, 10 bias.

## Prerequisites (once, per HPC account)

1. **Environment + package** (see the top-level pipeline README for detail):
   ```bash
   conda activate image_processing
   pip install -e /path/to/photometric_pipeline --no-deps --no-build-isolation
   ```
2. **Astrometry.net index files** — one shared read-only copy, pointed to via an
   environment variable (no hard-coded path overrides this):
   ```bash
   export CASSA_ASTROMETRY_INDEX=/path/to/astrometry_data
   ```
3. **Network** — the Phase 3 cross-match reaches APASS/Pan-STARRS/SDSS, so the
   compute node must have outbound internet for notebook `04`.

## Building the lecture PDF

```bash
cd slides
pdflatex cassa_photometry_lecture.tex   # run twice for the outline/TOC
```

Requires a TeX install with `beamer` (e.g. `texlive-latex-recommended`). Uses only
stock Beamer themes so it compiles anywhere.

## Running the hands-on notebooks

Start Jupyter in the `image_processing` env and run the notebooks **in order**.
Notebook `00` is a setup check and `01` is self-contained (it fabricates its own
example file); `02`–`05` each depend on the previous notebook's outputs.

Every notebook ends with **exercises that have worked solutions** — the point is
to read and re-run them, then change something and see what moves.

| Notebook | What it does |
|----------|--------------|
| `00_setup` | Verify env, package, index dir, `solve-field`, dataset. |
| `01_data_model_sci_err_dq` | Inventory the raw frames; build a toy MEF and decode SCI / ERR / DQ. |
| `02_phase1_calibration` | Run ISR → `calibrated_*.fits`; watch ERR/DQ populate. |
| `03_phase2_integration` | Align + stack + WCS-solve → `Master_*.fits`. |
| `04_phase3_photometry` | Zero point (online catalog) → catalog CSV + fluxcal. |
| `05_phase4_diagnostics` | Generate and read the diagnostics report. |

**Before running:** the config cell at the top of each notebook defaults to paths
*relative to the notebooks directory* — `../raw`, `../work`, and
`../../astrometry_data` — which work as-is if you launch Jupyter from
`workshop/notebooks/` with the layout above. On a shared machine where the index
files live elsewhere, edit that cell (the block marked `EDIT THESE PATHS`) or set
`CASSA_ASTROMETRY_INDEX` in your environment. The cell fails immediately with a
clear message if `RAW_DIR` does not exist.

> Notebooks 02–05 assume the provided dataset is a single target/filter for a clean
> walkthrough. For multi-filter datasets, repeat Phase 3 per filter.
