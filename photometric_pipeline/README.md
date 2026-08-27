# cassa-photometry

An end-to-end photometric reduction pipeline for the **CASSA Observatory /
iTelescope** network. It takes raw FITS frames and produces WCS-solved,
flux-calibrated images and source catalogs — with a **full error budget carried
from the raw pixels to the final magnitudes**.

The pipeline runs in three phases, all sharing one set of conventions
(configuration, logging, and multi-extension FITS I/O):

| Phase | Command | Input → Output |
|-------|---------|----------------|
| 1. Calibration (ISR) | `cassa-calibrate` | raw FITS → `calibrated_*.fits` (electrons) |
| 2. Integration | `cassa-integrate` | calibrated frames → WCS-solved `Master_*.fits` |
| 3. Photometry | `cassa-photometry` | master frames → `*_fluxcal.fits`, `*_catalog.csv`, `*_segmap.fits` |
| 4. Diagnostics | `cassa-diagnose` | all stage outputs → `diagnostics/diagnostics_report.pdf`, PNGs, `metrics.json` |

Plus `cassa-verify` (compare catalog magnitudes against APASS/Pan-STARRS/SDSS)
and `cassa-run` (all three phases chained).

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

```bash
# Phase 1: calibrate raw frames
cassa-calibrate -i /data/raw -o /data/calibrated

# Phase 2: align, stack, solve WCS (writes into a RunNN_* subdirectory)
cassa-integrate /data/calibrated --yes

# Phase 3: zero point, flux calibration, catalogs
cassa-photometry /data/calibrated/Run01_.../

# Verify against reference catalogs
cassa-verify /data/calibrated/Run01_.../

# Phase 4: diagnostics across all stages
cassa-diagnose /data/calibrated/Run01_.../ --raw /data/raw

# Or run everything end to end
cassa-run -i /data/raw -o /data/calibrated
```

## Phase 4 diagnostics

`cassa-diagnose` validates that the pipeline is behaving, stage by stage, and
writes a multi-page `diagnostics/diagnostics_report.pdf` plus per-stage PNGs, a
`metrics.json`, and a pipeline-health summary. Given a Phase 2 run directory it
auto-discovers each stage's product (and pairs them by filter); `--raw` points at
the raw frame/dir for stage 0. Each stage estimates the **PSF (FWHM via 2D
Gaussian fits + a stacked radial profile)** and adds stage-specific checks:

- **stage_0 (raw):** image + histogram, background, star detection, saturation map.
- **stage_1 (calibrated):** SCI/ERR/DQ panels, DQ flag counts, ERR-vs-signal (Poisson) check, background-flatness improvement vs raw.
- **stage_2 (master):** SCI/ERR/DQ, SNR map, depth boost vs a single frame (≈√N), WCS status + field centre/plate scale.
- **stage_3 (photometry):** Mag vs Mag_Error, SNR vs Mag, number counts / limiting magnitude, source map, star/galaxy shape split, ZP/`MAGZERR`.

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
```

```bash
cassa-photometry /data/masters --config my_config.yaml
```

## Development

```bash
pip install -e ".[dev]"
pytest
```
