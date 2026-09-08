# cassa-photometry

An end-to-end photometric reduction pipeline for the **CASSA Observatory**,
extensible to any imaging setup through an instrument profile. It takes raw FITS
frames and produces WCS-solved, flux-calibrated images and source catalogs —
with a **full error budget carried from the raw pixels to the final
magnitudes**.

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
  sigma-clipped), and writes `FLUXERR_ISO`, `MAGERR_ISO`, `SNR` and `FLAGS`
  columns in the catalog.

### Which magnitude to use

The catalog carries several, because they measure different things:

| Column | What it is | Use it for |
|--------|-----------|------------|
| **`MAG_BEST`** | PSF magnitude for point sources, Kron for extended | **the default choice** |
| `MAG_PSF` | PSF-fitted; the matched filter for a star | faint stars, crowded fields |
| `MAG_AUTO` | Kron elliptical (SExtractor's `AUTO`) | galaxies, total flux |
| `MAG_APER` | Measured in the zero point's own aperture | anything needing exact consistency with `MAGZERO` |
| `MAG_ISO` | **Isophotal** — flux above the detection threshold | comparison with older catalogs only |

`MAG_ISO` is not a total magnitude and never was. It captures a
brightness-dependent fraction of a source, so applying an aperture-derived zero
point to it produces a *tilt* rather than an offset. Measured against simulated
data with known truth:

```
              median error   scatter   trend with brightness
MAG_BEST         +0.034       0.091      +0.004 mag/mag
MAG_APER         -0.006       0.110      -0.011 mag/mag
MAG_ISO          +0.531       0.546      +0.348 mag/mag     <- a tilt
```

`MAG_ISO`'s error runs from +0.18 mag at V=12 to +1.70 mag at V=16.5. It is kept
under its own name for compatibility; the pipeline says so in the log.

Sources are classified as `STAR`, `EXTENDED`, `AMBIGUOUS`, `SATURATED` or
`EDGE` from `MAG_PSF - MAG_AUTO` against the frame's own stellar locus, with a
`CLASSLIM` recording the magnitude below which the classes stop being
meaningful. `AMBIGUOUS` is the honest answer, not a failure.

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
./install.sh
```

That is the whole procedure. `install.sh` picks an environment (an active conda
env, else conda if you have it, else a plain `python3 -m venv`), installs every
dependency **including a working plate solver**, installs the pipeline, and then
runs `cassa-doctor` to prove it worked. Running it again updates in place.

The conda route is the one to prefer: `environment.yml` pins Python, brings the
scientific stack as prebuilt binaries, and includes JupyterLab for the workshop
notebooks. Install [Miniforge](https://conda-forge.org/download/) first if you
have no conda.

On **Windows**, install [WSL](https://learn.microsoft.com/windows/wsl/install)
and run the same commands inside it. There is no native Windows install and
`install.ps1` does not attempt one: no plate solver is published for Windows —
conda-forge builds `astrometry` for `linux-64` and `osx-64` only, and the PyPI
in-process solver ships no Windows wheel — so a native environment would build
and then fail at the first WCS solve. Running `install.ps1` prints the WSL
steps.

```powershell
.\install.ps1        # prints the WSL instructions; installs nothing
```

### If something is wrong

```bash
cassa-doctor
```

It reports your Python and package versions against what the pipeline requires,
which plate-solving backends are usable, where astrometry index files will come
from, whether the reference catalogs are reachable, and whether the working
directory is writable — one line each, with the fix for anything that failed.
Paste its output when asking for help.

### About the plate solver

Phase 2 needs one of two backends, and the installer arranges whichever suits
your machine:

| Backend | Where it comes from | Available on |
|---|---|---|
| `solve-field` | conda-forge `astrometry`, or your system package manager | Linux x86-64, macOS Intel. Preferred when present. |
| in-process | `pip install "cassa-photometry[solver]"` | Linux x86-64, macOS Intel **and Apple Silicon**. |

Neither exists for Windows or for ARM Linux, which is why those go through WSL
and x86-64 respectively.

The solver is deliberately **not** in `environment.yml`: conda-forge has no
`osx-arm64` build of `astrometry`, and an environment file has no way to say
"only on some platforms", so listing it there makes `conda env create` fail
outright on an Apple Silicon Mac. `install.sh` installs it afterwards, picking
the backend the platform can run.

The two must never both be installed: conda's Astrometry.net package ships
Python bindings that import under the same name as the PyPI solver. `install.sh`
adds the PyPI one only when `solve-field` is absent.

### Astrometry index files

Plate solving matches star patterns against sky *index files*. **You do not need
to download them.** The pipeline works out which files a field requires and
fetches those from the public Astrometry.net server, caching them under
`~/.cache/cassa-photometry/astrometry`.

The saving is the point: a CASSA 8-inch field needs **4 files, 165 MB** out of
the ~34 GB the server offers — 0.5% of the set.

To see what a dataset will need before committing to the download, or to prepare
a laptop before going somewhere without a network:

```bash
cassa-index-fetch --from-headers raw/ --dry-run    # what it needs, and how big
cassa-index-fetch --from-headers raw/              # fetch it
```

If you already have a full local set, point the pipeline at it and nothing will
be downloaded:

```bash
export CASSA_ASTROMETRY_INDEX=/path/to/astrometry_data     # environment variable
# or set phase2.astrometry_index_dir in a config YAML
# otherwise the pipeline defaults to ./astrometry_data
```

### Development

```bash
./install.sh          # includes pytest and ruff
pytest
make lint
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

The *raw* tree is the one layout the pipeline does not own, so it reads whatever
the acquisition software wrote. A night arrives sorted by frame type and target:

```
raw/
  20260903/
    BIAS/untargeted/   20260903T054529.715412_untargeted_Blue_BIAS_....fits
    DARK/untargeted/   ...
    FLAT/untargeted/   ...
    LIGHT/m22/         20260903T043707.171323_m22_Luminance_LIGHT_....fits
```

Point `-i` at `raw/` (or at the whole archive, or at a single flat directory of
frames) and it is searched recursively. The folder names are documentation, not
classification: every frame is sorted by its `IMAGETYP`, `FILTER` and `EXPTIME`
header, so a frame filed under the wrong folder is still reduced as what it is.

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
- **stage_3 (photometry):** Mag vs MAGERR_ISO, SNR vs Mag, number counts / limiting magnitude, source map, star/galaxy shape split, ZP/`MAGZERR`.

## Verification (`cassa-verify`)

An independent check of the zero point: catalog magnitudes are cross-matched
against **APASS → Pan-STARRS → SDSS** and reported side by side, with **both**
error bars — your `MAGERR_ISO` (from the ERR plane) and the reference catalog's
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
| `cassa8` | CASSA 8-inch f/5 Newtonian + QHY miniCAM8M (IMX585 mono), LRGB+SHO wheel. |

Only the setups this observatory operates are built in, so the list is a
statement about CASSA rather than a directory of everyone's telescopes.
**Describing your own setup never requires editing this package** — see
[Adding a setup](#adding-a-setup).

Select one on any command, or in the config file:

```bash
cassa-calibrate -i /data/raw -o work/phase1 --instrument cassa8
cassa-run -i /data/raw -o work --instrument cassa8
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

Three routes, in increasing effort. Use the first that works.

**1. A `detector:` block — no code.** If your camera simply does not record its
own gain and read noise, say so and reduce with `generic`:

```yaml
# my_config.yaml
instrument: generic
detector:
  gain: 1.4                 # e-/ADU
  read_noise: 7.0           # e-
  saturation_adu: 60000
  pixel_scale_arcsec: 0.40  # unbinned
  filter_map: {Sloan-R: R_Photo}
  science_bands: {Sloan-R: R}
```

These fill in **where the header is silent**; the frame remains the first
authority on its own data. `override_header: true` inverts that, and says so in
the log every time it discards a card.

**2. A profile class in a local file.** For anything header-dependent — a CMOS
conversion-gain curve, a multi-amplifier layout. Copy
[`examples/profiles/template.py`](examples/profiles/template.py), then:

```yaml
instrument_module: my_profile.py:MyObservatoryProfile
```

No installation, no edit to this repository — which is what keeps your clone
rebasable against upstream.

**3. An installed package.** To distribute a profile, register it under the
`cassa_photometry.instruments` entry-point group in your own `pyproject.toml`.

`instruments/cassa.py` is the worked example in the package;
[`examples/profiles/`](examples/profiles/) holds a commented template and a
retired real-world profile.

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
`MAG_INST` with `MAG_ISO` as `NaN`.

LRGB `Red`/`Green`/`Blue` map to the R/G/B catalog bands and luminance to V.
These are imaging filters, not Johnson-Cousins or Sloan, so a colour term
remains — fine for differential photometry, worth stating for absolute work.

## Customizing a run

Three levels, cheapest first, all documented in
[`docs/CUSTOMIZING.md`](docs/CUSTOMIZING.md):

1. **Config** — turn any reduction step off, reorder the steps, add your own,
   change any threshold. What ran is recorded in the product (`CALPLAN` /
   `CALSKIP` in phase 1, `STEPPLAN` / `STEPSKIP` in phases 2–3,
   `step_plans` in phase 4's `metrics.json`), so a customised reduction cannot
   pass for a default one.
   ```yaml
   phase1:
     steps:
       cosmic_rays: false            # or: exclude: [cosmic_rays, flat]
       # reject cosmic rays BEFORE flat fielding rather than after
       order: [linearity, overscan, bad_pixel_mask, bias, dark,
               cosmic_rays, flat, measure_fwhm]
   phase3:
     steps:
       psf_photometry: false
       custom: {write_bright_list: local/my_steps.py:write_bright_list}
   ```
   An order that breaks a real dependency (flat before bias) is **refused when
   the config is read**, naming the violation — but genuine choices, like
   cosmic rays before or after the flat, are permitted. Every command takes
   `--skip`, `--only` and `--show-plan`; `cassa-run` also takes `--from`/`--to`.
   ```bash
   cassa-run -i raw -o work --show-plan     # print the plan, reduce nothing
   cassa-run -i raw -o work --from 2 --to 3 # resume without recalibrating
   ```
2. **Describe your setup** — a `detector:` block or your own instrument profile,
   without editing this package.
3. **Change an algorithm** — documented extension seams, and how to keep a
   modified clone mergeable with upstream.

## Working offline

The pipeline is meant to be usable from a laptop away from the university.

```bash
cassa-index-fetch --from-headers raw/    # astrometry indexes, once
cassa-run -i raw -o work --instrument cassa8   # populates the catalog cache
# ...later, with no network:
cassa-photometry work/phase2 --offline
```

Reference-catalog queries are cached under `~/.cache/cassa-photometry/catalogs`,
so a field reduced once can be re-reduced with no network at all. `--offline`
never touches the network and says clearly when a field is not in the cache,
rather than hanging on a series of timeouts.

## Simulated data with known truth

```bash
cassa-simulate --preset workshop --out sim
cassa-run -i sim/raw -o sim/work --instrument cassa8
```

Generates a complete observing run — bias, dark, flat and science frames — from
a real star field (Gaia DR3 positions and magnitudes, so the frames genuinely
plate-solve), and writes every source's true magnitude and every frame's true
seeing, transparency and zero point alongside. The dataset is distributed as a
**seed rather than a download**: the same command gives everyone the same
frames.

This is how the numbers quoted in this README were measured, and it is the
honest way to check that a change to the pipeline improved something rather than
merely altered it.

## Configuration

All tunable parameters (cosmic-ray thresholds, stacking sigma-clip, aperture
radii, zero-point match tolerance, saturation level, …) live in
`cassa_photometry.config`. Override any subset with a YAML file:

```yaml
# my_config.yaml
instrument: cassa8               # see "Instrument profiles" above
log_level: INFO

phase1:
  saturation_adu: 60000          # fallback only; see "Bad pixels" below
  fallback_gain: 1.0             # last resort; used only if header AND profile
  fallback_read_noise: 10.0      # are both silent, and warned about when it is
  allow_precalibrated: false     # reduce frames with a non-empty CALSTAT anyway
  bpm_dark_rate_factor: 20.0     # hot if dark current > 20x the median rate
  bpm_bias_noise_factor: 5.0     # unstable if bias scatter > 5x the median
  cr_max_fraction: 0.01          # discard a cosmic-ray mask claiming more than
                                 # 1% of a frame; no cosmic-ray rate reaches it,
                                 # so it is a crowded field being flagged as one
phase2:
  solver: auto                   # auto | solve-field | astrometry-py
  epoch_bin: night               # none | night | "6h" -- see below
  stack_weight: point_source     # point_source | extended
  fwhm_reject_factor: 1.6        # drop frames softer than 1.6x the median
phase3:
  zp_sigma_clip: 2.5
  offline: false                 # reference catalogs from cache only
  aperture_r_factor: 2.0         # aperture radius in units of the MEASURED FWHM
  match_radius_sigma: 3.0        # cross-match radius = sigma * ASTRMS ...
  match_radius_min_arcsec: 1.0   # ... clamped to this floor ...
  match_radius_max_arcsec: 5.0   # ... and this cap
```

Note that `phase3.fwhm` is now only a **fallback**. Apertures are sized from
`FWHMPX`, which phases 1 and 2 measure from the stars in the frame itself; a
fixed aperture whatever the seeing was is worth up to 0.3 mag of per-epoch error
in the zero point.

### Epochs

Frames are grouped by observing night, binned on **local noon** rather than the
UT date — an observing night crosses UT midnight at most longitudes, and binning
on the date string cuts one night into two half-depth masters.

```yaml
phase2:
  epoch_bin: night     # default
  epoch_bin: none      # all dates together (the pre-0.2 behaviour)
  epoch_bin: "6h"      # sub-night bins, for fast variables
```

Nightly binning averages intra-night variability. That is right for a supernova
and wrong for a short-period variable, which wants a sub-night bin.

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
