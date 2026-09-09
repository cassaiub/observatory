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

> **Looking something up?** [`docs/REFERENCE.md`](docs/REFERENCE.md) has every
> command and flag, every configuration key with its default, every environment
> variable, every FITS keyword and catalog column, and the Python API.
> [`docs/CUSTOMIZING.md`](docs/CUSTOMIZING.md) covers changing what a run does.

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

The WCS carries an uncertainty too, and **every backend produces it** — the
pipeline measures the residual itself rather than taking the solver's word,
because `solve-field` reports matched star pairs and ASTAP does not. Where a
backend supplies its own matched stars they are used; otherwise the frame's
detected sources are projected through the solved WCS and paired against the
reference catalog Phase 3 already caches.

| Keyword | Meaning |
|---------|---------|
| `CRDER1`, `CRDER2` | FITS-standard random error per axis, in degrees |
| `ASTRMS` | total RMS residual, in arcsec |
| `ASTNSTAR` | number of stars matched |
| `ASTRMSRC` | which route produced it — the solver's own table, or the reference catalog |

This distinguishes a solve that *converged* from one that actually **fits**, and
it is what `cassa-verify` uses to size its cross-match radius. Phase 4 reports it
in the stage 2 panel and in `metrics.json`.

Measuring it in-pipeline is what keeps the two install paths equivalent. If the
residual came from the solver, a frame solved by ASTAP would carry none, Phase 3
would silently fall back to a fixed 2″ match radius, and the photometry would
depend on which solver a given machine happened to install.

## Installation

```bash
git clone https://github.com/cassaiub/observatory.git
cd observatory/photometric_pipeline
./install.sh
```

That is the whole procedure on **Linux and macOS**. `install.sh` picks an
environment (an active conda env, else conda if you have it, else a plain
`python3 -m venv`), installs every dependency **including a working plate
solver**, installs the pipeline, and then runs `cassa-doctor` to prove it
worked. Running it again updates in place.

The conda route is the one to prefer: `environment.yml` pins Python, brings the
scientific stack as prebuilt binaries, and includes JupyterLab for the workshop
notebooks. Install [Miniforge](https://conda-forge.org/download/) first if you
have no conda.

### Platform support

| Platform | Route | Solver you get |
|---|---|---|
| **Linux x86-64** | `./install.sh` | ASTAP, else `solve-field`, else in-process |
| **Linux aarch64** (Raspberry Pi, ARM servers) | `./install.sh --conda` | **ASTAP only** — nothing else is published for ARM Linux |
| **macOS Intel** | `./install.sh` | ASTAP, else `solve-field`, else in-process |
| **macOS Apple Silicon** | `./install.sh` | ASTAP, else in-process (`solve-field` has no `osx-arm64` build) |
| **Windows** | WSL, then the Linux route | as Linux x86-64 |

### Step by step, per operating system

**Linux** (Debian/Ubuntu shown; any distribution works)

```bash
sudo apt update && sudo apt install -y git curl unzip
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh          # accept the defaults
exec $SHELL -l                                    # pick up conda on PATH

git clone https://github.com/cassaiub/observatory.git
cd observatory/photometric_pipeline
./install.sh
conda activate cassa-photometry
cassa-doctor
```

On **aarch64** pass `--conda`: the in-process solver publishes no ARM Linux
wheel, so the conda route plus ASTAP is the combination that works.

**macOS** (Intel and Apple Silicon alike — the installer name is built from your
own machine, so the commands are identical)

```bash
xcode-select --install                            # once, for git and the toolchain
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
exec $SHELL -l

git clone https://github.com/cassaiub/observatory.git
cd observatory/photometric_pipeline
./install.sh
conda activate cassa-photometry
cassa-doctor
```

Gatekeeper quarantines an unsigned ASTAP binary downloaded from the web. If
`cassa-doctor` reports the binary as present but unusable, clear the attribute:

```bash
xattr -d com.apple.quarantine "$(command -v astap_cli)"
```

**Windows** — through WSL, which is real x86-64 Linux. From an **administrator**
PowerShell:

```powershell
wsl --install                    # reboot when prompted
```

Then open **Ubuntu** from the Start menu, finish creating your user, and follow
the Linux instructions above inside that window. Every command is the Linux one,
and `cassa-doctor` reports native Windows as a failure with this same guidance
rather than letting a run die later at the solver.

```powershell
.\install.ps1        # prints the WSL instructions; installs nothing
```

Keep the repository and the work directory in your WSL home. Your Windows drives
are mounted under `/mnt/c`, so data downloaded on the Windows side is reachable,
but reducing a night across `/mnt` is several times slower.

> **Why WSL and not a native install?** ASTAP does publish a Windows build, so
> the solver is no longer the obstacle it once was. What is missing is the rest:
> `install.sh` is a shell script, the test suite and the packaging classifiers
> target Linux and macOS, and no native Windows configuration has been verified
> end to end. WSL is the route that is tested, so it is the route that is
> documented.

### If something is wrong

```bash
cassa-doctor
```

It reports your platform, your Python and package versions against what the
pipeline requires, which plate-solving backends are usable, how many ASTAP star
tiles are cached, where astrometry index files will come from, whether the
reference catalogs are reachable, and whether the working directory is writable
— one line each, with the fix for anything that failed. Paste its output when
asking for help.

### About the plate solver

Phase 2 needs one, and **three backends** can provide it. `install.sh` tries
them in this order and the first that installs wins; `solver: auto` then prefers
whichever is present, in the same order.

| Backend | Where it comes from | Available on |
|---|---|---|
| **`astap`** | `apt install astap-cli`, else the upstream zip | Linux **x86-64 and aarch64**, macOS **Intel and Apple Silicon** |
| `solve-field` | conda-forge `astrometry`, or your system package manager | Linux x86-64, macOS Intel |
| `astrometry-py` | `pip install "cassa-photometry[solver]"` | Linux x86-64, macOS Intel and Apple Silicon |

ASTAP leads because it is the only backend that exists everywhere this pipeline
runs — on ARM Linux there is no astrometry.net solver from either channel — and
because it is an 875 KB binary with no Python dependency. On the CASSA test
frames it solved in **0.11 s** against 41 s for the in-process solver, and it
recovered the true plate scale from a header that stated it wrongly.

`solve-field` and the PyPI solver must **never both** be installed: conda-forge's
Astrometry.net package ships Python bindings that import under the same name as
the PyPI solver, and whichever lands second wins. `install.sh` adds the PyPI one
only when `solve-field` is absent; `cassa-doctor` reports which backend is in
use.

The solver is deliberately **not** in `environment.yml`: conda-forge has no
`osx-arm64` build of `astrometry`, and an environment file has no way to say
"only on some platforms", so listing it there makes `conda env create` fail
outright on an Apple Silicon Mac. `install.sh` installs it afterwards, picking
the backend the platform can run.

**The products do not depend on which one ran.** See
[Astrometric error](#astrometric-error) above.

### Star databases and index files: nothing to download

Plate solving matches star patterns against pre-computed sky data. **You do not
need to fetch any of it in advance.** Whichever backend is in use, the pipeline
reads the frame — pointing from the header, field size from the plate scale and
the array dimensions — works out which files that field requires, checks what is
already cached, and fetches only the difference.

| | ASTAP | Astrometry.net |
|---|---|---|
| Per field | **~6 MB** of star tiles | ~246 MB of index files |
| Full published set | 859 MB | ~34 GB |
| Cache | `~/.cache/cassa-photometry/astap` | `~/.cache/cassa-photometry/astrometry` |
| Config override | `phase2.astap_db_dir` | `phase2.index_cache_dir` |
| Environment variable | `CASSA_ASTAP_DB` | `CASSA_INDEX_CACHE` |
| Turn fetching off | `phase2.astap_db_download: false` | `phase2.index_download: false` |

ASTAP publishes no per-tile URL — its star databases are distributed as one
859 MB ZIP. The pipeline reads individual tiles out of that archive with HTTP
range requests: a ZIP's central directory lists every member's offset, so only
the members a field needs are ever transferred. The central directory costs
85 KB, the tiles about 6 MB, and the observatory hosts nothing.

Repointing therefore costs only the tiles the new field adds. Measured on the
workshop data: 6.3 MB and 65 s for the first field from an empty cache, then
**1.2 s** for the next field with no download at all.

For the Astrometry.net path you can see or pre-fetch the cost:

```bash
cassa-index-fetch --from-headers raw/ --dry-run    # what it needs, and how big
cassa-index-fetch --from-headers raw/              # fetch it
```

A full local set you manage yourself always wins over the cache, and nothing is
downloaded:

```bash
export CASSA_ASTROMETRY_INDEX=/path/to/astrometry_data   # a full index set
export CASSA_ASTAP_DB=/path/to/astap_database            # an installed ASTAP db
```

Setting `phase2.astap_db_download: false` (or `index_download: false`) still
performs the selection and reports exactly which files are missing, rather than
hanging — the air-gapped case.

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
# Reduce the field once, on a network. This populates every cache:
cassa-run -i raw -o work --instrument cassa8
# ...later, with no network at all:
cassa-photometry work/phase2 --offline
```

Three caches make this work, all under `~/.cache/cassa-photometry/`: the solver's
sky data (`astap/` or `astrometry/`) and the reference-catalog queries
(`catalogs/`). A field reduced once can be re-reduced with no network.

To prepare a laptop for a field it has never seen, warm the solver cache in
advance:

```bash
cassa-index-fetch --from-headers raw/ --dry-run   # astrometry.net: the cost first
cassa-index-fetch --from-headers raw/             # ...then fetch it
```

ASTAP has no equivalent pre-fetch command because it does not need one: its
tiles are ~6 MB per field, fetched during the solve itself. Run the reduction
once while you have a connection and the cache is warm.

`--offline` never touches the network and says clearly when a field is not in the
cache, rather than hanging on a series of timeouts. `phase2.astap_db_download:
false` and `phase2.index_download: false` do the same for the solver: selection
still runs, and the run names exactly which files are missing.

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
  solver: auto                   # auto | astap | solve-field | astrometry-py
  astap_db_series: auto          # auto | d50 | d05 | g05 -- chosen by field size
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

## Documentation

| Document | What it covers |
|---|---|
| This README | What the pipeline is, and how to run it. |
| [`docs/REFERENCE.md`](docs/REFERENCE.md) | Every command, flag, config key, environment variable, FITS keyword, catalog column and public function. |
| [`docs/CUSTOMIZING.md`](docs/CUSTOMIZING.md) | Excluding, reordering and adding steps; describing your setup; changing an algorithm and staying mergeable. |
| [`docs/CHANGELOG.md`](docs/CHANGELOG.md) | What changed, and what it means for existing data. |
| [Observatory manual](../docs/main.pdf) | The full treatment: the algorithms, the data-flow diagrams, and the other two observatory packages. |
| [`workshop/`](workshop/) | A lecture, a participant handbook, and five notebooks that reduce a night with known truth. |

## Development

```bash
./install.sh          # includes pytest and ruff
pytest                # 474 tests
make lint
```

The test suite runs without a network and without a plate solver installed: the
subprocess boundary is faked where that is what is under test, and anything
needing outbound access is marked `network` (`pytest -m "not network"`).
