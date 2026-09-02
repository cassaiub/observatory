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

## Instrument profiles

Everything the pipeline needs to know about a telescope lives behind one object,
the **instrument profile**: how the detector reports its gain and read noise, how
its filters are named, which frames are calibrations. Every phase resolves its
profile by name, so the same pipeline reduces an amateur rig and a professional
observatory without a fork.

| Name | Setup |
|------|-------|
| `generic` | **Default.** Standard FITS keywords only — correct for any setup whose acquisition software writes them (see the header spec). |
| `itelescope` | The iTelescope hosted network (T11, T24, T32, T68), with per-telescope detector constants and its luminance-flat-for-red convention. |
| `cassa8` | CASSA 8-inch f/5 Newtonian + QHY miniCAM8M (IMX585 mono), LRGB+SHO wheel. |

Select one on any command, or in the config file:

```bash
cassa-calibrate -i /data/raw -o work/phase1 --instrument cassa8
cassa-run -i /data/raw -o work --instrument itelescope
```

```yaml
# my_config.yaml
instrument: cassa8
```

`--instrument` overrides the config file, which overrides the `generic` default.

### What a profile is for

The header is always the first authority: a frame carrying `EGAIN` and
`READNOIS` is believed by every profile. A profile exists to supply what the
header *cannot* say — and to refuse to guess when it does not know:

1. the frame's own header, always;
2. the profile's hardware knowledge, when the header is silent;
3. `None` — reported as a warning naming the frame, then the configured
   `phase1.fallback_gain` / `fallback_read_noise`, so an assumed error budget is
   never mistaken for a measured one.

### Adding a setup

Subclass `InstrumentProfile`, override only what differs, and add one line to
`instruments/registry.py`. The base class is concrete, so a setup that writes
standard keywords may need to override nothing at all beyond its detector
constants. `instruments/cassa.py` is a worked example.

### CMOS cameras: gain is not a number

On a CMOS detector the system gain in e-/ADU and the read noise are *functions of
the gain setting*, not constants — the IMX585 in the CASSA 8-inch has an HCG/LCG
transition near gain 30 where read noise drops from ~3.2 e- to ~1.3 e-. A single
scalar cannot describe this, so `Cassa8InchProfile` interpolates a curve keyed on
the `GAIN` card. `get_gain(header)` receives the whole header, so this needs no
special interface.

The same profile deliberately does **not** read `GAIN` as e-/ADU, because on this
camera `GAIN` is the unitless setting — reading it as e-/ADU is a ~100x error.
Write `EGAIN` at acquisition and the question never arises.

The mode comes from `READOUTM` when the header carries it. Without it, the mode is
*inferred* from the gain setting against `HCG_SWITCH_GAIN` — right for a camera
left on auto, a guess for one driven manually. Writing `READOUTM` removes the
guess:

```
GAIN=29, no READOUTM      -> LCG   0.1500 e-/ADU   3.20 e-
GAIN=29, READOUTM=HCG     -> HCG   0.1400 e-/ADU   1.30 e-   <- same gain card
```

> **The `cassa8` detector curves are provisional** — published IMX585 figures, not
> measurements of this camera. Measure it with the sibling `camera_characterization`
> package, whose PTC stage produces exactly this curve, then replace
> `DETECTOR_CURVES` — `Cassa8InchProfile.set_detector_curve(curve, mode="LCG")`
> loads one mode at a time — and set `MEASURED = True`.

### Subframes, binning, and already-calibrated frames

Two geometry situations reach phase 1, and they need opposite answers.

A **windowed camera** produces science frames smaller than the full-frame master
calibrations. That is a routine ROI setup, so the masters — bias, dark, flat and
the bad-pixel mask — are cropped to the frame's region using the
`XORGSUBF`/`YORGSUBF` origin the profile reports, uncertainty planes included.

A **binning mismatch** cannot be cropped away, and the frame is skipped. The two
used to produce the same misleading "Check binning" message; now each says which
it is.

Phase 1 also refuses frames that report having been calibrated already — a
non-empty `CALSTAT`, or data whose `BUNIT` says electrons rather than ADU.
Reducing a reduced frame yields a plausible-looking image with a wrong error
budget and nothing downstream can detect it, so the check runs once per file
during the header scan:

```
[WARNING] Skipping 1 frame(s) that report prior calibration: sci_done.fits (CALSTAT='BDF')
[WARNING] Raw frames should have an empty CALSTAT and be in ADU. Set
          phase1.allow_precalibrated: true to reduce them anyway.
```

### Filters without a reference band

Narrowband filters (Ha, SII, OIII) and blocked wheel positions have no broadband
counterpart in APASS, Pan-STARRS or SDSS. Cross-matching one anyway yields a zero
point that is numerically valid and physically meaningless, so `science_band()`
returns `None` for them and phase 3 **skips the zero point and flux calibration**.
Detection still runs and the catalog is still written, carrying
`Instrumental_Mag` with `Absolute_Mag` as `NaN`.

LRGB `Red`/`Green`/`Blue` map to the R/G/B catalog bands and luminance to V.
These are imaging filters, not Johnson-Cousins or Sloan, so a colour term
remains — fine for differential photometry, worth stating for absolute work.

## Configuration

All tunable parameters (cosmic-ray thresholds, stacking sigma-clip, aperture
radii, zero-point match tolerance, saturation level, …) live in
`cassa_photometry.config`. Override any subset with a YAML file:

```yaml
# my_config.yaml
instrument: cassa8               # see "Instrument profiles" above
phase1:
  saturation_adu: 60000          # fallback only; see "Bad pixels" below
  fallback_gain: 1.0             # last resort; used only if header AND profile
  fallback_read_noise: 10.0      # are both silent, and warned about when it is
  allow_precalibrated: false     # reduce frames with a non-empty CALSTAT anyway
  bpm_dark_rate_factor: 20.0     # hot if dark current > 20x the median rate
  bpm_bias_noise_factor: 5.0     # unstable if bias scatter > 5x the median
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

### Bad pixels

The bad-pixel mask is built from every master available, because each
calibration product sees a defect class the others are blind to:

| Master | Catches | Why the others miss it |
|--------|---------|------------------------|
| flats  | dead, low-QE, dust | — |
| darks  | hot pixels | a few seconds of dark current is a rounding error against the lamp signal in a flat |
| biases | unstable / flickering | these read out at the right *level*, just never the same one twice |

Thresholds are relative to each master's own median rather than absolute in ADU
or e⁻/s, so the defaults carry from one detector to the next without retuning.
Set a `*_factor` to `0` to disable that test; a master that is absent is simply
skipped. Missing hot pixels is not harmless: astroscrappy re-detects them as
cosmic rays in *every* frame, which is how a one-off event ends up flagged in
all of them.

Saturation is a property of the detector, so `Phase1Config.saturation_adu` is
only the fallback. A profile that knows its camera's full well should override
`InstrumentProfile.get_saturation` (the base implementation already reads
`SATURATE` / `SATLEVEL` / `FULLWELL` from the header), so that adding a
telescope does not mean shipping a config file with it.

### Flat combination order

Each flat is divided by its own median **before** the stack is combined, not
after. Twilight flats fade as they are taken — on the workshop set the level
drifts 27–57% within a single filter — so a median across raw ADU frames tracks
the fading sky rather than the pixel response, and the stack scatter that
becomes the master's `ERR` measures that drift too.

Combining first inflates the master flat's relative uncertainty from ~0.11% (the
photon-noise floor) to ~4.3%. Since `flat_correct` propagates *relative* error,
that lands hardest on bright pixels: error bars on 5000–20000 e⁻ sources were
inflated by up to 3.5×, which feeds straight into Phase 3 magnitude errors and
zero-point weighting. The flat field itself also differs by up to 2.5%.

### Dark subtraction and exposure times

Each master dark records the exposure it represents (plus gain, combine count,
and whether the bias pedestal was removed). Science frames of a *matching*
exposure get the dark subtracted unscaled — the preferred case, since amp glow
is non-linear and does not survive scaling cleanly. A mismatch is rescaled by
the exposure ratio and logged.

Refusing to scale does not enforce a matched exposure; it only hides the
mismatch. A 120s master dark subtracted whole from a 60s frame removes twice the
thermal signal the frame accumulated, uniformly and silently. Scaling needs a
bias-subtracted dark (a raw dark is `bias + rate·t`, and scaling that scales the
pedestal too), so the pipeline warns if asked to rescale one that still has its
pedestal. Darks of mixed exposure times in one directory are normalised to a
common exposure before combining rather than averaged incoherently.

## Development

```bash
pip install -e ".[dev]"
pytest
```
