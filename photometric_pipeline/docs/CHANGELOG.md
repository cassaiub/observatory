# Changelog

Notable changes to `cassa-photometry`. Entries that change a **delivered
number** are marked **[science]** and say what moved and by how much — a
photometric pipeline that changes its answers silently is not usable for
research.

Products carry a `CALVERS` card recording the calibration vintage they were
written with (`fits_utils.CALVERS`). Files with no `CALVERS` are vintage 1.
Catalogs of different vintages must not be combined without re-reduction.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — 0.2.0.dev0

Work tracked in `docs/master-plan.md`. This release is a scientific-hardening
pass across all four phases, a self-provisioning installation, and an
instrument layer built around the CASSA 8-inch.

### Added

- **ASTAP as a third plate-solving backend, and the one the installer prefers.**
  `install.sh` now tries **ASTAP → `solve-field` → in-process** and uses the
  first that installs; `phase2.solver: auto` prefers them in the same order.
  ASTAP is an 875 KB binary with no Python dependency, packaged as `astap-cli`
  on Debian and Ubuntu and published for Linux x86-64 and aarch64 and for macOS
  on both Intel and Apple Silicon — it is the only backend that exists on ARM
  Linux, where neither conda-forge nor PyPI ships an astrometry.net solver at
  all. On the CASSA test frames it solved in **0.11 s** against 41 s for the
  in-process solver, and it recovered the true plate scale from a header that
  stated it wrongly. New settings: `astap_path`, `astap_db_dir`,
  `astap_db_series`, `astap_db_url`, `astap_db_download`,
  `astap_tile_neighbours`. `WCSSOLVR` records which backend solved a frame.
- **`ASTRMS` is now measured by the pipeline, so every backend produces it.**
  ASTAP reports no matched star pairs, and phase 3 sizes its reference-catalog
  cross-match radius from `ASTRMS` — so taking the residual from the solver
  would have made the photometry depend on which solver a machine happened to
  install (a missing `ASTRMS` silently falls back to a fixed 2″ radius). The new
  `phase2_integration.astrometry_qc` measures it from things every backend
  produces: the solved WCS, sources detected with `sep`, and reference
  positions. Where a backend supplies its own catalog stars they are used (no
  network); otherwise the cached reference catalog is reused. The header records
  which route was taken as **`ASTRMSRC`**. Measured end to end: ASTAP solved,
  then the pipeline measured `ASTRMS 0.358″` from 83 Gaia stars.
- **ASTAP star tiles are fetched per field, like the astrometry indexes.** ASTAP
  publishes no per-tile URL — its databases are one 859 MB ZIP — so
  `astap_db.py` reads individual tiles out of that archive with HTTP range
  requests: a ZIP's central directory lists every member's offset, so only the
  members a field needs are transferred. Measured: central directory 85 KB, then
  **6.2 MB in 33 range requests** for the NGC 7331 field, against 859 MB for the
  whole archive. Repointing costs only the tiles the new field adds — first
  field 65 s from a cold cache, **1.2 s** for the next with no download.
  `astap_db_download: false` still performs the selection and names the missing
  tiles, for an air-gapped machine. Nothing is hosted by the observatory.
- **A native Windows installer, provisionally.** Windows was unsupported
  because no plate solver was published for it. ASTAP ended that: it ships
  command-line builds for `win64`, `win32` and `win11_aarch64`, and all 16
  runtime dependencies have Windows wheels or are pure Python. `install.ps1`
  now mirrors `install.sh` — it picks an environment, installs the package,
  fetches the ASTAP build matching the machine's architecture, and runs
  `cassa-doctor` — instead of printing WSL instructions and exiting.
  `install.ps1 -Wsl` still prints those. The packaging declares
  `Operating System :: Microsoft :: Windows`, and `cassa-doctor` reports the
  platform as **WARN, not FAIL** — native Windows has not been verified end to
  end, so failing would assert a problem nobody has observed while passing
  would assert a guarantee nobody has earned. WSL remains the recommended route
  until someone confirms it; `docs/WINDOWS-TESTING.md` is the checklist.

  One trap found while building this, worth recording because it is invisible:
  SourceForge serves a *browser* an HTML "your download will start shortly"
  page instead of the file, and PowerShell's `Invoke-WebRequest` identifies as
  a browser by default — so the naive download lands 114 KB of HTML named
  `astap.zip`. The installer asks with a non-browser user agent (322 KB,
  beginning `PK`) and verifies the ZIP magic before trusting it.
- **A complete reference document, `docs/REFERENCE.md`**, and Part VI of the
  observatory manual: every command and flag, every configuration key with its
  default, every environment variable, every FITS keyword read and written,
  every catalog column, the DQ bits, and the module-by-module Python API.
- **[science] A cosmic-ray mask that cannot be real is now thrown away.** In a
  crowded field astroscrappy cannot separate a stellar peak from a cosmic ray,
  and it *replaces* what it flags: on a globular cluster it deleted 23% of the
  flux at the flagged pixels, 3.2% of all the light in the frame, and the DQ
  flags it left behind then disqualified those same stars as zero-point
  calibrators in phase 3. Phase 1 now measures the share of the frame each mask
  claims and discards it above `phase1.cr_max_fraction` (default 1%), recording
  `CRVETO` in the header. The threshold is physical, not tuned: cosmic rays
  arrive at a few events per cm^2 per minute, well under 0.1% of the pixels even
  for a large sensor in a long exposure. Reject cosmic rays in the phase 2 stack
  instead, where they do not repeat between frames.
- **The raw tree from the acquisition software is read as delivered.** A night
  arrives sorted into `<date>/BIAS|DARK|FLAT|LIGHT/<target>/`, and every command
  that takes a raw directory (`cassa-calibrate -i`, `cassa-run -i`,
  `cassa-diagnose --raw`, `cassa-integrate --raw-dir`, `cassa-index-fetch
  --from-headers`) now searches it recursively, through the one helper
  `paths.find_raw_frames`. A flat directory of frames — what `cassa-simulate`
  writes — reads exactly as before. Classification is unchanged and still comes
  from `IMAGETYP`/`FILTER`/`EXPTIME`, never from a folder name, so a misfiled
  frame is reduced as what it is. Two raw frames sharing a filename in different
  directories no longer overwrite each other's output: the colliding names carry
  their raw subdirectory, and the run says so.
- **[science] The catalog now carries magnitudes that are actually total.**
  `MAG_APER` (the zero point's own aperture), `MAG_AUTO`/`KRON_RADIUS` (Kron
  elliptical), `MAG_PSF`/`CHI2_PSF` (PSF-fitted) and `MAG_BEST` (PSF for point
  sources, Kron for extended). Measured against simulated truth, for stars
  brighter than 16.5:

  | column | median error | scatter | trend with brightness |
  |---|---|---|---|
  | `MAG_BEST` | +0.034 mag | 0.091 | +0.004 mag/mag |
  | `MAG_APER` | −0.006 mag | 0.110 | −0.011 mag/mag |
  | `MAG_ISO` (unchanged) | +0.531 mag | 0.546 | **+0.348 mag/mag** |

  `MAG_ISO` keeps its meaning and its name — it is isophotal, and its error runs
  from +0.18 mag at V=12 to +1.70 at V=16.5. **Use `MAG_BEST` or `MAG_APER`.**
- **[science] Star/galaxy classification from concentration, not shape.**
  `CLASS`, `CLASS_STAR`, `PSF_MINUS_AUTO`, `FLUX_RADIUS`, `FWHM_IMAGE`, with the
  threshold calibrated per frame from the stellar locus and a reported
  `CLASSLIM` below which no class is asserted. The previous
  `ELLIPTICITY < 0.15` test measured shape: a face-on elliptical galaxy passed
  as a star, a trailed star did not — and a galaxy in the zero point biases
  every magnitude in the frame. The zero point and `cassa-verify` now both cut
  on `CLASS`.
- **[science] An aperture correction, measured per frame from a curve of
  growth.** `MAGZERO` is now a *total-flux* zero point. The correction is fitted
  as a spatial **surface** where enough stars support it (94–178 per frame on
  CASSA data), because an 8-inch Newtonian's PSF broadens off-axis and a single
  number is right at the centre and wrong in the corners; where it cannot be, the
  measured spread is folded into `MAGZERR` rather than discarded. Recorded as
  `APCOR`, `APCORMTH`, `APCORRMS`.
- **[science] Gaia DR3 synthetic photometry as the primary reference catalog**,
  ahead of APASS→Pan-STARRS→SDSS. It supplies Johnson-Cousins *and* SDSS
  magnitudes for the same stars, which fixes a defect neither plan had caught:
  `FILTER='R'` is Johnson R and was being calibrated against Sloan r′, a
  **−0.21 ± 0.036 mag** systematic measured on real stars in the NGC 7331 field.
  ~487 calibrators per field, and the product header records which system was
  used.
- **Reference-catalog caching and `--offline`.** Every successful query is
  persisted, so a field reduced once can be re-reduced with no network at all.
- **`targets.yaml` and header-declared targets.** `TARGNAME`/`OBJTYPE`/
  `TARGRA`/`TARGDEC`/`HOSTGAL`/`DISCDATE` route a frame; a `--targets` file
  overrides and supplements them. A frame declaring nothing is a field, exactly
  as before. A frame declaring a supernova is warned about loudly, because host
  removal is not implemented yet and an ordinary measurement of one is biased low.
- **Reduction steps can be switched off from config.** Each phase has a
  `steps:` block; a skipped step is recorded in the product (`CALSKIP`) and in
  the log, so a partially-reduced frame cannot pass for a finished one.
- **`docs/CUSTOMIZING.md`** — the four levels of customization (config,
  instrument, code, package), the documented extension seams, and how to keep a
  modified clone mergeable with upstream.
- **CI** (`.github/workflows/ci.yml`) — lint plus the suite with **no network**
  on Python 3.10/3.11/3.12, and a separate job that runs `./install.sh` from a
  clean checkout and reduces simulated data with it. "Clone and run one command"
  is a claim, so it is exercised rather than asserted.
- **`cassa_photometry.psf`** — PSF measurement moved to the top level (all three
  science phases depend on it now), with `build_epsf` for an empirical PSF from
  a frame's own stars. The old import path still works.
- **Astrometry index files fetch themselves.** A field needs a handful of the
  ~34 GB the Astrometry.net server offers, and the pipeline now works out which
  and downloads those — from the public server, so nothing has to be hosted and a
  fresh clone needs no configuration. Measured: a CASSA 8-inch field needs **4
  files, 165 MB, 0.5% of the set**. A populated `astrometry_index_dir` still
  wins, so existing installations are unchanged.
  - `cassa-index-fetch` prefetches for a pointing or for every frame in a
    directory (`--from-headers`), with `--dry-run` to see the cost first.
  - Selection samples the search cone's *boundary*, not just its centre: 143 of
    400 random pointings straddle a tile edge, and centre-only selection silently
    misses the neighbour they need.
- **Two plate-solving backends behind one interface.** `solve-field` when it is
  installed, and an in-process solver (the PyPI `astrometry` package) otherwise,
  so `pip install` alone gives a working Phase 2. `phase2.solver` forces one;
  `cassa-doctor` reports which is in use.
- **`cassa-simulate` — simulated CASSA 8-inch data with known truth.** Generates
  a complete run (bias, dark, flat, science) through a forward model that is the
  deliberate inverse of phase 1, and writes the true magnitude of every source
  and the true seeing, transparency and zero point of every frame alongside. The
  dataset is distributed as a **seed, not a download**: `cassa-simulate --preset
  workshop` gives everyone the same frames.
  - The star field is **real** — Gaia DR3 positions with synthetic
    Johnson-Cousins magnitudes — so `solve-field` genuinely solves the frames
    and the reference-catalog cross-match genuinely finds the stars. Verified: a
    workshop-sized frame solves with 67 matches of 90 index stars, 0 conflicts,
    recovering the field centre to 0.13".
  - The PSF is a **Moffat that broadens off-axis**, not a Gaussian. A Gaussian
    has an aperture correction of essentially zero, which is why the existing
    synthetic tests passed while the aperture correction was missing.
- **Three ways to describe an observing setup, none of which touch this
  package.** A `detector:` config block (gain, read noise, saturation, pixel
  scale, filter tables) for the common case of a camera that does not record its
  own constants; `instrument_module: my_profile.py:MyProfile` to load a profile
  class from a local file; and a `cassa_photometry.instruments` entry-point group
  for a profile shipped by another package. `examples/profiles/` holds a
  commented template and a retired real-world profile, both covered by tests.
- **`InstrumentProfile.get_pixel_scale`** — resolves the plate scale from the
  header (`SECPIX` and friends), else pixel size and focal length, else the
  profile's optics with binning applied, else an existing WCS. Phase 2 needs it
  for a solver scale hint and phase 3 for its catalog search cone; the workshop
  frames carry `SECPIX = 0.4` that nothing read, so `PIXSCALE` came out `0.0`.
- **`InstrumentProfile.apply_linearity`** — a hook for detector non-linearity,
  identity by default. No profile measures one yet.
- **One-command installation.** `./install.sh` (Linux/macOS, venv or conda)
  builds the environment, installs the package, ensures a working plate solver,
  and verifies the result. `environment.yml` is now a real file in the repo
  rather than copy-paste text in the manual.
  *(Corrected 2026-09-08: this entry originally described `install.ps1` as
  building a conda environment on Windows. It never did, and the script now
  installs nothing, printing the WSL setup steps instead. Linux and macOS are
  the supported platforms; see the 2026-09-09 note under Added — ASTAP has since
  made a Windows solver available, but the install tooling and test coverage
  have not followed, so WSL remains the documented route.)*
- **`cassa-doctor`** — reports Python, package versions against their declared
  floors, which plate-solver backends are usable, the astrometry index
  situation, network reachability, and write permissions. Non-zero exit only on
  a genuine failure, so it works at the end of an install and in CI.
- **`cassa`** — an umbrella command (`cassa calibrate`, `cassa run`, …)
  dispatching to the same functions the individual `cassa-*` scripts use.
- `Phase4Config` — phase 4 had no config block at all.
- `CALVERS` calibration-vintage card, stamped on every product
  (`fits_utils.CALVERS`, currently **2**).
- `docs/master-plan.md` — the authoritative, resumable plan for this release.
- `Makefile` (`install`, `test`, `lint`, `doctor`, `docs`, `clean`) and a ruff
  configuration; `network` and `slow` pytest markers.

### Changed

- Version bumped to `0.2.0.dev0`.
- **Config files are now validated.** A wrong type (`phase3: {fwhm: "big"}`) or
  a mistyped key is reported against its full dotted path instead of being
  applied silently and failing later somewhere unrelated.
- **The order of authority between header, config and profile is now
  explicit.** `InstrumentProfile` separates `gain_from_header` / 
  `read_noise_from_header` / `saturation_from_header` — what the frame itself
  claims — from `get_gain` / `get_read_noise` / `get_saturation` — what the
  profile concludes. Without that seam a user-configured value was shadowed by
  the profile's own fallback table, because a profile like `cassa8` never
  returns `None`. The order is now header → `detector:` → profile → fallbacks.
- **`log_level` now works.** It was declared in the config and read by nothing;
  every logger was hardcoded to `INFO`.
- **Phases 1, 3 and 4 write a run log**, as phase 2 always has. A batch run
  previously left no record of what it did.
- **[science] Instrumental magnitudes are normalised by exposure time**, using
  the true mean single-frame exposure (`EXPMEAN`). A master is a weighted *mean*,
  so dividing by `TOT_EXP` would be wrong by a factor of N, and the group key's
  5-second quantisation carries a 0.02 mag systematic of its own.
- **[science] Flux-calibrated images are `Jy/pixel`, not `Jy`**, and apply the
  band's AB−Vega offset. The conversion is the AB relation but B/V/R/I are
  calibrated against Vega magnitudes: without the offset B fluxes were ~8.6%
  wrong. Recorded as `PHOTSYS`, `ABOFFSET`, `ZPAB`, `PHOTWAVE`.
- **[science] The reference-catalog search cone covers the whole field.** It used
  the half-*width* and `CD1_1`; on a real master that asked for 5.04′ where 7.13′
  was needed, leaving half the field area with no calibrators.
- **[science] Zero-point calibrators are quality-cut.** Stars with a saturated,
  bad or cosmic-ray pixel in the aperture are rejected (now possible at all,
  since phase 2 propagates DQ); the background annulus is a sigma-clipped
  *median* moved outside the PSF wings; matching is mutually-nearest and uses the
  radius derived from `ASTRMS` rather than a hardcoded 2″.
- **[science] Non-detections are kept.** Sources with non-positive flux stay in
  the catalog with a finite uncertainty and a 3σ `LIMIT_MAG`, instead of being
  dropped — which biased faint number counts and made upper limits impossible.
- **A zero-point failure no longer destroys the catalog.** They shared one
  `try`/`except`, so a network blip left an empty phase3 directory.
- **[science] Plate solving: three defects fixed, each of which was costing
  solves.**
  - *The frame's own pointing now outranks everything.* The solver was handed
    the **previous group's** solved centre in preference to the frame's own
    `OBJCTRA`/`OBJCTDEC`, so any run covering more than one field pointed at the
    wrong sky for every group after the first.
  - *A failed solve now widens its hints instead of discarding them*, then drops
    the scale hint entirely as a last resort. A **wrong** scale makes a solve
    fail where a missing one only makes it slow — and the workshop frames claim
    `SECPIX = 0.4` where blind solving measures **0.591"/px**, an error of 48%.
  - *The extractor's noise level is measured, not assumed.* `--sigma` is the
    assumed image noise in ADU, not a threshold multiplier, so the hardcoded `8`
    was specific to one telescope's data. On CASSA-like frames (true noise 56
    ADU) it returned **499 sources of which 28 were real**.
- **Phase 4 pairs its products by filter.** The last-resort catalog lookup took
  the alphabetically first file, so a report could describe a B-band catalog
  beside a V-band master without saying so. The PDF is also closed even when a
  stage raises, rather than being left truncated.
- **The solved plate scale is now checked against the header's claim** and a
  disagreement is reported. Free, and it is how a wrong `SECPIX` gets noticed
  instead of quietly breaking every solve on that instrument.
- **Phase 2 no longer blocks on `input()`.** `HardwareManager.allocate_resources`
  prompted from inside `setup()`, *before* the `--yes` check, so any run without
  a terminal — `cassa-run`, a batch job, a notebook, a test — died with
  `EOFError`. It now prompts only when a terminal is attached, and its output
  goes to the log rather than `print`.
- `--version` works on every command; only `cassa-calibrate` had it.
- An unknown `--instrument` is reported by the command the user typed, rather
  than surfacing as a `KeyError` from inside a phase.
- Dependency floors set where an API break is real (`astropy>=6.0`,
  `photutils>=2.0`, `astroscrappy>=1.1`, `ccdproc>=2.4`); `requests` and `sep`
  added. `astroscrappy` mattered: versions ≤1.0 returned the cleaned array
  already in electrons, which would silently square the gain correction.
- `docs/astrometry-index-plan.md` and `docs/target-driven-photometry-plan.md`
  are now superseded reference documents; progress is tracked in
  `docs/master-plan.md`.

### Verified

Against simulated data with known truth, end to end on a 2048×1400 CASSA 8-inch
field (`cassa-simulate --preset workshop`):

- All three filters plate-solve on the first pass using **4 index files
  (165 MB)** of the ~34 GB available, at 0.41–0.51″ astrometric RMS over 59–74
  matched stars, recovering the plate scale to **0.03%**.
- `MAG_BEST` recovers truth to **+0.034 mag** with **no trend** against
  brightness (+0.004 mag/mag).
- **131 of 144** true stars classified `STAR`; the galaxy `EXTENDED`.
- Phase 1 recovers the injected dark current exactly (0.0050 e⁻/s) and the
  measured `FWHMPX` tracks the true delivered PSF to **5%**.
- An `--offline` re-run reproduces all three zero points **exactly** from cache.
- `./install.sh --venv` builds a working installation from a clean checkout;
  the suite passes there (311 tests) as well as under conda (313).

On real iTelescope data (NGC 7331, B/V/R, 1024², 0.17° field), for the solver
work:

- ASTAP and the in-process solver recover plate scales agreeing to **0.089%**
  (0.5905 vs 0.5910″/px), from a header that states 0.4″/px — a 32% error, which
  both the ASTAP warning and `solve_scale_warn_frac` flagged.
- A cold ASTAP cache fetched **10 tiles, 6.3 MB, 33 range requests**, solved in
  0.11 s, and the pipeline then measured `ASTRMS 0.358″` from 83 Gaia stars.
- A second field with the cache warm: **1.2 s, no download**.
- The suite is now **474 tests** (was 439), including 27 covering the tile
  arithmetic, the cache and the ASTAP backend, and 8 covering the residual — the
  cos(dec) de-projection, mutual pairing, and the refusal to report a residual
  from a single star.

### Fixed

- **The in-process solver was completely broken, and so was its residual.** Two
  pre-existing bugs, both found while wiring up ASTAP, both fixed:
  1. `astrometry` 4.3.0 declares `Solver(index_files: list[pathlib.Path])` and
     calls `path.resolve()`; `inprocess.py` passed `str()`, so every solve raised
     `AttributeError: 'str' object has no attribute 'resolve'`. With
     `astrometry>=4.3` pinned in `pyproject.toml`, **the pip-only path — the
     thing that makes `pip install` sufficient — failed for every user.**
  2. `_matched()` read `star.metadata["x"]` from `Match.stars`, but those are
     *catalog* stars whose metadata comes from the index file and never
     contained our pixel coordinates — so `residuals()` always returned `None`
     and `ASTRMS` was never written. It now projects the extracted sources
     through the solved WCS and pairs them with the catalog by position.

  Both were masked by two `test_solvers.py` failures that stopped at "no backend
  available". Verified: the solve now completes in 41.9 s and the residual
  returns `(0.801, 0.664, 45)` where it previously returned `None`.

### Removed

- **The `itelescope` profile.** The pipeline now ships only the setups CASSA
  operates. This changes no existing number: the profile's hardware table is
  keyed on `TELESCOP`, and the workshop frames report `CDK700`, which matched no
  entry — so it had always fallen through to its `DEFAULT` of gain 1.0 /
  read noise 10.0, identical to `generic` with the configured fallbacks. The
  class is kept, annotated, as `examples/profiles/itelescope.py` and still
  loadable via `instrument_module:`.

## [0.1.0] — 2026-08

Initial release: four-phase reduction (ISR → integration/astrometry →
photometry/zero point → diagnostics) with SCI/ERR/DQ error propagation, the
pluggable instrument layer, and the CASSA 8-inch profile.
