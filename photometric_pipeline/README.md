# cassa-photometry

An end-to-end photometric reduction pipeline for the **CASSA Observatory /
iTelescope** network. It takes raw FITS frames and produces WCS-solved,
flux-calibrated images and source catalogs — with a **full error budget carried
from the raw pixels to the final magnitudes**.

The pipeline runs in four phases, all sharing one set of conventions
(configuration, logging, and multi-extension FITS I/O). Each phase reads the
previous phase's directory and writes its own:

| Phase | Command | Input → Output | Directory |
|-------|---------|----------------|-----------|
| 1. Calibration (ISR) | `cassa-calibrate` | raw FITS → `calibrated_*.fits` (electrons) | `work/phase1` |
| 2. Integration | `cassa-integrate` | calibrated frames → WCS-solved `Master_*.fits` | `work/phase2` |
| 3. Photometry | `cassa-photometry` | master frames → `*_fluxcal.fits`, `*_catalog.csv`, `*_segmap.fits` | `work/phase3` |
| 4. Diagnostics | `cassa-diagnose` | all stage outputs → `diagnostics_report.pdf`, PNGs, `metrics.json` | `work/phase4` |

A re-run replaces that phase's products in place rather than accumulating copies.

Plus `cassa-verify` (compare catalog magnitudes against APASS/Pan-STARRS/SDSS)
and `cassa-run` (phases 1–3 chained).

## Error propagation & data model

Every science image is a **multi-extension FITS** file with three planes:

```
[0] PRIMARY  SCI  — science data (electrons; Jy after flux calibration)
[1] ERR           — 1-sigma per-pixel uncertainty (same units as SCI)
[2] DQ            — integer bitmask (saturated / bad pixel / cosmic ray / no-data)
```

- **Phase 1** seeds the error budget with `ccdproc.create_deviation` (Poisson +
  read noise) and builds master calibrations *with* uncertainty, so `ccdproc`
  propagates the ERR plane through bias/dark/flat/gain automatically.
- **Phase 2** warps and scales the ERR plane alongside the data and materializes
  the inverse-variance-weighted stack variance as the master ERR plane.
- **Phase 3** uses the ERR plane for aperture and segment photometry, computes a
  **filter-wise zero point with uncertainty** (`MAGZERO` / `MAGZERR`,
  sigma-clipped), and writes `Flux_Error`, `Mag_Error`, `SNR` and `FLAGS`
  columns in the catalog.

### Astrometric error

The WCS carries an uncertainty too. After a successful solve, Phase 2 reads
`solve-field`'s `.corr` table — the detected stars paired with their
index-catalog counterparts — and records the RMS of those residuals:

| Keyword | Meaning |
|---------|---------|
| `CRDER1`, `CRDER2` | FITS-standard random error per axis, in degrees |
| `ASTRMS` | total RMS residual, in arcsec |
| `ASTNSTAR` | number of stars matched against the index |

This distinguishes a solve that *converged* from one that actually **fits**, and
it is what `cassa-verify` uses to size its cross-match radius. Phase 4 reports it
in the stage 2 panel and in `metrics.json`.

## Installation

```bash
git clone <repo-url>
cd cassa_observatory/photometric_pipeline
pip install .            # or: pip install -e .   (editable, for development)
```

This installs all Python dependencies (astropy, ccdproc, photutils, astroalign,
astroquery, …) and the `cassa-*` commands.

### System dependency: Astrometry.net (`solve-field`)

Phase 2 shells out to `solve-field`, which is **not** a Python package and cannot
be installed with `pip`. Install it separately:

```bash
# conda (recommended, cross-platform)
conda install -c conda-forge astrometry

# Debian/Ubuntu
sudo apt-get install astrometry.net
```

### Astrometry index files

`solve-field` needs sky index files (several GB). They are **not** bundled and are
git-ignored. Download the appropriate `index-*.fits` files for your field of
view and point the pipeline at them via any of:

```bash
export CASSA_ASTROMETRY_INDEX=/path/to/astrometry_data     # environment variable
# or set phase2.astrometry_index_dir in a config YAML
# otherwise the pipeline defaults to ./astrometry_data
```

## Usage

Each phase writes into its own directory under one work directory:

```
work/
  phase1/   calibrated frames          (cassa-calibrate)
  phase2/   master stacks + WCS        (cassa-integrate)
  phase3/   flux-calibrated + catalogs (cassa-photometry)
  phase4/   diagnostics report         (cassa-diagnose)
```

A re-run replaces that phase's products in place rather than accumulating copies.

```bash
# Phase 1: calibrate raw frames
cassa-calibrate -i /data/raw -o /data/work/phase1

# Phase 2: align, stack, solve WCS (writes to /data/work/phase2)
cassa-integrate /data/work/phase1 --yes

# Phase 3: zero point, flux calibration, catalogs (writes to /data/work/phase3)
cassa-photometry /data/work/phase2

# Verify against reference catalogs
cassa-verify /data/work/phase3

# Phase 4: diagnostics across all stages (writes to /data/work/phase4)
cassa-diagnose /data/work/phase2 --raw /data/raw

# Or run everything end to end (-o is the work directory)
cassa-run -i /data/raw -o /data/work
```

Each phase also accepts an explicit `-o/--output` (or `--outdir`) if you want a
different location.

## Phase 4 diagnostics

`cassa-diagnose` validates that the pipeline is behaving, stage by stage, and
writes a multi-page `diagnostics_report.pdf` plus per-stage PNGs, a
`metrics.json`, and a pipeline-health summary into `work/phase4`. Given the
Phase 2 directory it auto-discovers each stage's product from the sibling phase
directories (and pairs them by filter); `--raw` points at the raw frame/dir for
stage 0. Each stage estimates the **PSF (FWHM via 2D
Gaussian fits + a stacked radial profile)** and adds stage-specific checks:

- **stage_0 (raw):** image + histogram, background, star detection, saturation map.
- **stage_1 (calibrated):** SCI/ERR/DQ panels, DQ flag counts, ERR-vs-signal (Poisson) check, background-flatness improvement vs raw.
- **stage_2 (master):** SCI/ERR/DQ, SNR map, depth boost vs a single frame (≈√N), WCS status + field centre/plate scale, astrometric RMS (`ASTRMS`) and matched-star count.
- **stage_3 (photometry):** Mag vs Mag_Error, SNR vs Mag, number counts / limiting magnitude, source map, star/galaxy shape split, ZP/`MAGZERR`.

## Verification (`cassa-verify`)

An independent check of the zero point: catalog magnitudes are cross-matched
against **APASS → Pan-STARRS → SDSS** and reported side by side, with **both**
error bars — your `Mag_Error` (from the ERR plane) and the reference catalog's
own uncertainty — so a disagreement can be judged against the errors rather than
eyeballed:

```
Cross-match radius: 2.09"  [3 x ASTRMS 0.696" (WCS fit to 26 stars)]

--- APASS VERIFICATION REPORT (Filter: r'mag | match radius: 2.09") ---
Obj ID  | RA (deg)   | DEC (deg)  | Your Mag  | Your Err  | APASS Mag  | APASS Err  | Delta
172     | 339.19497  | 34.46303   | 11.969    | 0.027     | 12.010     | 0.083      | -0.041 mag
```

The match radius is **derived from the astrometric solution** rather than fixed:
`match_radius_sigma × ASTRMS`, clamped to `[match_radius_min_arcsec,
match_radius_max_arcsec]`. The floor matters — reference positions carry their
own error and stars have moved since the catalog epoch, so even an excellent
solve should not use a radius of a few tenths of an arcsecond. When a master
carries no `ASTRMS` (solved before this was recorded, or a failed solve), it
falls back to the fixed `phase3.zp_match_tol_arcsec` and says so.

A catalog that reports no uncertainty for a star shows `n/a` rather than an
invented number — Pan-STARRS does this for stars saturated in *its* survey.

## Configuration

All tunable parameters (cosmic-ray thresholds, stacking sigma-clip, aperture
radii, zero-point match tolerance, saturation level, …) live in
`cassa_photometry.config`. Override any subset with a YAML file:

```yaml
# my_config.yaml
phase1:
  saturation_adu: 60000
phase3:
  fwhm: 4.0
  zp_sigma_clip: 2.5
  match_radius_sigma: 3.0        # cassa-verify radius = sigma * ASTRMS ...
  match_radius_min_arcsec: 1.0   # ... clamped to this floor ...
  match_radius_max_arcsec: 5.0   # ... and this cap
```

```bash
cassa-photometry /data/masters --config my_config.yaml
```

## Development

```bash
pip install -e ".[dev]"
pytest
```
