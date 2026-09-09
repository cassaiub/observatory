# Reference

Everything the pipeline exposes, in one place: every command and flag, every
configuration key with its default, every environment variable, every FITS
keyword it reads and writes, every catalog column, and the Python API behind
them.

This is the lookup document. For *why* the pipeline does what it does, read the
[observatory manual](../../docs/main.pdf) (Part III); for how to change what it
does, read [`CUSTOMIZING.md`](CUSTOMIZING.md).

- [1. Commands](#1-commands)
- [2. Configuration](#2-configuration)
- [3. Environment variables](#3-environment-variables)
- [4. The step registry](#4-the-step-registry)
- [5. FITS keywords](#5-fits-keywords)
- [6. Catalog columns](#6-catalog-columns)
- [7. Data quality (DQ) bits](#7-data-quality-dq-bits)
- [8. Python API](#8-python-api)
- [9. What `cassa-doctor` checks](#9-what-cassa-doctor-checks)

---

## 1. Commands

Ten console entry points. `cassa <command>` is an umbrella for all of them, so
`cassa calibrate` and `cassa-calibrate` are the same program.

| Command | Phase | Reads | Writes |
|---|---|---|---|
| `cassa-calibrate` | 1 | raw FITS tree | `work/phase1/calibrated_*.fits` |
| `cassa-integrate` | 2 | `work/phase1` | `work/phase2/Master_*.fits`, QA PDF |
| `cassa-photometry` | 3 | `work/phase2` | `work/phase3/*_fluxcal.fits`, `*_catalog.csv`, `*_segmap.fits` |
| `cassa-diagnose` | 4 | all phases | `work/phase4/diagnostics_report.pdf`, PNGs, `metrics.json` |
| `cassa-run` | 1→3 | raw FITS tree | `work/phase1`, `phase2`, `phase3` |
| `cassa-verify` | — | `*_catalog.csv` | a report on stdout |
| `cassa-doctor` | — | your machine | a report on stdout |
| `cassa-index-fetch` | — | frame headers | the astrometry index cache |
| `cassa-simulate` | — | a seed | a complete simulated run with truth files |
| `cassa` | — | — | dispatches to any of the above |

### Options shared by every phase command

| Flag | Effect |
|---|---|
| `-c`, `--config FILE` | YAML file overriding any subset of the defaults. |
| `--instrument NAME` | Instrument profile. Overrides the config file. One of `generic`, `cassa8`, or anything registered by a plugin. |
| `--version` | Print the version and exit. |
| `-h`, `--help` | Print usage and exit. |

### Options shared by phases 1–4

| Flag | Effect |
|---|---|
| `--skip STEP` | Skip a reduction step. Repeatable. |
| `--only STEP` | Run **only** these steps, skipping every other one in the phase. Repeatable. |
| `--show-plan` | Print the resolved step order and exit, reducing nothing. |

### Resource options (`cassa-integrate`, `cassa-run`)

| Flag | Effect |
|---|---|
| `--cpu-fraction F` | Fraction of the machine's cores to use, `0`–`1`. |
| `--cores N` | Hard ceiling in cores, applied after the fraction. |
| `--max-memory-gb GB` | Memory the run should stay inside. Advisory — the estimate is checked against it and a run that will not fit says so before it starts. |

Any of these also **suppresses the interactive CPU prompt**: a stated budget is
an answer, and being asked again would make the flag useless in a batch job or a
notebook. They override the corresponding `phase2` config keys.

A step may be named bare (`--skip flat`) or qualified (`--skip phase1.flat`).
The qualified form only matters for `cassa-run`, where two phases both declare
`measure_fwhm`: there, a bare name that is ambiguous is an **error** rather than
a guess.

### `cassa-calibrate` — Phase 1, instrument signature removal

```
cassa-calibrate -i INPUT -o OUTPUT [-c CONFIG] [--instrument NAME]
                [--skip STEP] [--only STEP] [--show-plan]
```

| Flag | Effect |
|---|---|
| `-i`, `--input DIR` | **Required.** Raw FITS directory, searched recursively. The acquisition software's night tree (`<date>/BIAS\|DARK\|FLAT\|LIGHT/<target>/`) can be given as-is; frames are sorted by `IMAGETYP`/`FILTER`/`EXPTIME`, never by folder name. |
| `-o`, `--output DIR` | **Required.** Directory for calibrated frames. |

### `cassa-integrate` — Phase 2, alignment, stacking, WCS

```
cassa-integrate DATA_DIR [-o OUTPUT] [--keep-temps] [-y] [--raw-dir DIR]
                [-c CONFIG] [--instrument NAME] [--skip STEP] [--only STEP]
                [--show-plan]
```

| Flag | Effect |
|---|---|
| `DATA_DIR` | **Required.** Directory of calibrated frames (phase 1's output). |
| `-o`, `--output DIR` | Output directory. Default: the `phase2` directory beside the input. |
| `--keep-temps` | Keep the solver's temporary files instead of removing them. |
| `-y`, `--yes` | Skip the interactive confirmation of the resource estimate. |
| `--raw-dir DIR` | Raw frame directory for the QA PDF's raw panel. Only needed when the raws have moved since phase 1 stamped `RAWDIR`. |

### `cassa-photometry` — Phase 3, zero point, flux calibration, catalogs

```
cassa-photometry INPUT_PATH [--filter BAND] [--fwhm PX] [--threshold SIGMA]
                 [--outdir DIR] [--offline] [--targets FILE] [-c CONFIG]
                 [--instrument NAME] [--skip STEP] [--only STEP] [--show-plan]
```

| Flag | Effect |
|---|---|
| `INPUT_PATH` | **Required.** A master FITS file or a directory of them. |
| `--filter BAND` | Fallback science band (`R`/`G`/`B`/`V`/`I`) when the header cannot say. Default `R`. |
| `--fwhm PX` | FWHM estimate in pixels. A **fallback only** — apertures are sized from the measured `FWHMPX`. |
| `--threshold SIGMA` | Detection threshold in sigma. Default 5. |
| `--outdir DIR` | Output directory. Default: alongside the inputs. |
| `--offline` | Never query reference catalogs; use the local cache only, and say so clearly when a field is not in it. |
| `--targets FILE` | A `targets.yaml`. Frames carrying `TARGNAME`/`OBJTYPE` need no file. |

### `cassa-verify` — check the zero point against reference catalogs

```
cassa-verify INPUT_PATH [--filter BAND] [-c CONFIG] [--instrument NAME]
```

| Flag | Effect |
|---|---|
| `INPUT_PATH` | **Required.** A `_catalog.csv` file or a directory of them. |
| `--filter BAND` | Fallback reference-catalog column. Default `rmag`. |

### `cassa-diagnose` — Phase 4, per-stage diagnostics

```
cassa-diagnose [RUN_DIR] [--raw PATH] [--calibrated PATH] [--master PATH]
               [--catalog PATH] [--fluxcal PATH] [--outdir DIR] [-c CONFIG]
               [--instrument NAME] [--skip STEP] [--only STEP] [--show-plan]
```

| Flag | Effect |
|---|---|
| `RUN_DIR` | A phase 2 run directory. Everything else is auto-discovered from the sibling phase directories and paired by filter. |
| `--raw PATH` | Raw frame, or a raw directory searched recursively (stage 0). |
| `--calibrated PATH` | A `calibrated_*.fits` frame (stage 1). |
| `--master PATH` | A `Master_*.fits` stack (stage 2). |
| `--catalog PATH` | A `*_catalog.csv` (stage 3). |
| `--fluxcal PATH` | A `*_fluxcal.fits`, for the zero-point header. |
| `--outdir DIR` | Where the report goes. |

### `cassa-run` — phases 1 → 2 → 3, chained

```
cassa-run -i INPUT -o OUTPUT [--filter BAND] [--keep-temps] [--targets FILE]
          [-c CONFIG] [--from N] [--to N] [--instrument NAME]
          [--skip STEP] [--only STEP] [--show-plan]
```

| Flag | Effect |
|---|---|
| `-i`, `--input DIR` | **Required.** Raw FITS directory, searched recursively. |
| `-o`, `--output DIR` | **Required.** Work directory; phases write to `<output>/phase1`, `phase2`, `phase3`. |
| `--from {1,2,3}` | First phase to run. Default 1. Resumes a run whose later phases failed, without recalibrating. |
| `--to {1,2,3}` | Last phase to run. Default 3. |
| `--filter BAND` | Fallback science band for phase 3. |
| `--targets FILE` | A `targets.yaml` for phase 3. |

Phase 4 is deliberately not chained: it is a report on a finished run, not a
reduction step.

### `cassa-doctor` — check that this machine can run the pipeline

```
cassa-doctor [-c CONFIG] [--no-network]
```

| Flag | Effect |
|---|---|
| `--no-network` | Skip the reachability checks, for an air-gapped machine. |

Exit status is `1` if any check failed and `0` otherwise, so
`cassa-doctor && echo ok` works in a script. See
[section 9](#9-what-cassa-doctor-checks).

### `cassa-index-fetch` — pre-fetch astrometry index files

```
cassa-index-fetch [--ra DEG] [--dec DEG] [--scale ARCSEC] [--size PX]
                  [--from-headers DIR] [--dry-run] [--wide]
                  [-c CONFIG] [--instrument NAME]
```

| Flag | Effect |
|---|---|
| `--from-headers DIR` | Read pointing and scale from every FITS frame in a directory and fetch the union of what they need. |
| `--ra`, `--dec`, `--scale`, `--size` | Describe one field by hand instead. |
| `--dry-run` | List the files and the total size without downloading. |
| `--wide` | Include the finest (largest) indexes, which the pipeline otherwise fetches only if a first solve fails. |

This command is for the astrometry.net backend. **ASTAP needs no equivalent**:
its tiles are ~6 MB per field and are fetched during the solve itself.

### `cassa-simulate` — a complete dataset with known truth

```
cassa-simulate [--preset NAME] [-o OUT] [--seed N] [-c CONFIG]
```

| Flag | Effect |
|---|---|
| `--preset {workshop,tiny,multinight}` | Which dataset. Default `workshop`. |
| `-o`, `--out DIR` | Output directory. Default `./sim`. |
| `--seed N` | RNG seed. The dataset is reproducible from it, so everyone gets bit-identical frames. |

Writes `raw/` plus `truth_sources.csv`, `truth_frames.csv` and
`truth_config.yaml`, so every measurement can be **scored** rather than only
inspected.

---

## 2. Configuration

Every key below is valid in a YAML file passed with `-c`. Any subset may be
given; anything omitted keeps its default. A wrong type or a mistyped name is
reported against its full dotted path at load time rather than silently ignored.

### Top level

| Key | Default | Meaning |
|---|---|---|
| `instrument` | `generic` | Instrument profile name. `--instrument` overrides it. |
| `instrument_module` | `null` | A profile class in a local file, as `path/to/my_profile.py:MyProfile`. Takes precedence over `instrument`. |
| `detector` | *(empty)* | Detector constants for a setup with no profile — see below. |
| `targets_file` | `null` | A `targets.yaml` describing what is being observed. |
| `log_level` | `INFO` | `DEBUG`, `INFO`, `WARNING`, `ERROR`. |
| `phase1`…`phase4` | — | The per-phase blocks below. |

### `detector:` — describe a setup without writing a profile

Values fill in **where the header is silent**, preserving the order of
authority: header → this block → profile → configured fallbacks.

| Key | Default | Meaning |
|---|---|---|
| `gain` | `null` | e⁻/ADU. |
| `read_noise` | `null` | e⁻. |
| `saturation_adu` | `null` | Raw ADU at which pixels clip. |
| `pixel_scale_arcsec` | `null` | Unbinned arcsec/pixel. |
| `filter_map` | `{}` | Extra raw-filter → stacking-label entries, merged over the profile's map. |
| `science_bands` | `{}` | Extra raw-filter → science-band (`R`/`G`/`B`/`V`/`I`) entries. |
| `uncalibrated_filters` | `[]` | Filters to exclude from photometric calibration entirely. |
| `flat_proxies` | `{}` | filter → stand-in filter, for a setup that reuses one flat for another. |
| `override_header` | `false` | Let these values beat the frame header rather than fill in behind it. Logged loudly every time a card is discarded. |

### `phase1:` — instrument signature removal

| Key | Default | Meaning |
|---|---|---|
| `steps` | all on | See [section 4](#4-the-step-registry). |
| `master_sigma_clip` | `3.0` | Sigma-clip threshold when combining calibration frames. `0` disables rejection. |
| `apply_linearity` | `true` | Correct detector non-linearity in raw ADU. Only acts when the profile supplies a curve. *(Deprecated duplicate of `steps.linearity`.)* |
| `calibration_temp_tolerance_c` | `3.0` | How far a calibration frame's sensor temperature may differ from the science frames' before it is reported. |
| `cr_sigclip` | `4.5` | `astroscrappy.detect_cosmics` Laplacian threshold. |
| `cr_sigfrac` | `0.3` | Fractional threshold for neighbouring pixels. |
| `cr_objlim` | `5.0` | Contrast limit between the Laplacian and the fine-structure image. |
| `cr_max_fraction` | `0.01` | Largest share of a frame the cosmic-ray mask may claim before it is discarded as implausible. `1.0` accepts any mask. |
| `saturation_adu` | `50000.0` | Fallback saturation in raw ADU. A profile that knows its detector overrides it. |
| `fallback_gain` | `1.0` | Last-resort e⁻/ADU, used only when neither header nor profile can say. Warned about, by frame. |
| `fallback_read_noise` | `10.0` | Last-resort read noise in e⁻, same conditions. |
| `allow_precalibrated` | `false` | Reduce frames reporting prior calibration (non-empty `CALSTAT`, or data already in electrons). |
| `bpm_flat_low` | `0.5` | Flat sensitivity floor, as a fraction of the flat's own **local** level. |
| `bpm_flat_high` | `1.5` | Flat sensitivity ceiling, same basis. |
| `bpm_flat_smooth_px` | `51` | Box size of the median smoothing that defines "local". `0` thresholds against the normalised flat directly. |
| `bpm_dark_rate_factor` | `20.0` | Hot if the dark current exceeds this multiple of the median dark rate… |
| `bpm_dark_sigma` | `5.0` | …or this many robust sigmas above it, whichever cut is higher. |
| `bpm_bias_noise_factor` | `5.0` | Unstable if bias frame-to-frame scatter exceeds this multiple of the median scatter… |
| `bpm_bias_noise_sigma` | `5.0` | …or this many robust sigmas above it, whichever is higher. |

Set any `*_factor` to `0` to disable that test. A master that is absent is
simply skipped.

### `phase2:` — alignment, stacking, astrometry

**Registration and stacking**

| Key | Default | Meaning |
|---|---|---|
| `steps` | all on | See [section 4](#4-the-step-registry). |
| `align_detection_sigma` | `1.5` | Source-detection threshold for `astroalign`. |
| `align_min_area` | `4` | Minimum connected pixels for an alignment source. |
| `scale_min` / `scale_max` | `0.65` / `1.5` | Accepted frame-to-frame flux-scale range. Outside it, the frame is rejected. |
| `background_box` | `(50, 50)` | 2D background estimator box, in pixels. |
| `background_filter` | `(3, 3)` | Median filter applied to the background mesh. |
| `subtract_background` | `true` | Subtract a 2D background before stacking. Removes real extended emission with the sky, so extended-source work may want it off. The removed level is recorded as `BKGLEVEL` either way. *(Deprecated duplicate of `steps.subtract_background`.)* |
| `stack_sigma` | `3.0` | Sigma-clip threshold when combining. |
| `stack_maxiters` | `3` | Sigma-clip iterations. |
| `sigma_clip_min_frames` | `5` | Below this frame count nothing is rejected — too little information to identify an outlier. |
| `min_frames_to_stack` | `3` | Fewest aligned frames worth stacking at all. |
| `mask_bad_pixels` | `true` | Exclude pixels phase 1 flagged from the combine instead of averaging them in. |
| `stack_weight` | `point_source` | `point_source` uses 1/(σ·FWHM)², optimal for stars; `extended` uses 1/σ². |
| `anchor_sharpness_power` | `1.0` | Anchor ranking is SNR / FWHM^power. `0` ranks on SNR alone. |
| `fwhm_reject_factor` | `1.6` | Reject frames softer than this multiple of the group median. `0` keeps every frame. |
| `scale_aperture_factor` | `3.0` | Aperture for the flux-scale measurement, in units of the worse of the two frames' FWHM. Generous on purpose: a Moffat's wings mean a small aperture reads seeing as transparency. |
| `scale_aperture_fwhm_default` | `4.5` | Fallback FWHM when neither frame reports one. |
| `warp_order` | `3` | Spline order for the geometric warp. |
| `epoch_bin` | `night` | `none` (all dates together), `night` (binned on **local noon**, so a night crossing UT midnight stays whole), or a duration such as `"6h"` for sub-night bins. |
| `raw_dir` | `null` | Raw tree for the QA PDF's raw panel. Normally unset — phase 1 stamps `RAWDIR`. |
| `cpu_fraction` | `null` | Fraction of the machine's cores to use, `0`–`1`. `null` prompts on a terminal and takes 0.5 otherwise. Setting it suppresses the prompt. |
| `max_cores` | `null` | Hard ceiling in cores, applied after the fraction. For a share of a shared box rather than a share of the hardware. |
| `max_memory_gb` | `null` | Memory the run should stay inside, in GB. Advisory: nothing can cap what numpy allocates, but the estimate is checked against it and a run that will not fit says so up front. |

**Plate solving**

| Key | Default | Meaning |
|---|---|---|
| `solver` | `auto` | `auto`, `astap`, `solve-field`, `astrometry-py`. `auto` tries them in that order and uses the first that can run. |
| `astap_path` | `null` | Path to the ASTAP executable. `null` looks for `astap_cli`, then `astap`, on `PATH`. |
| `astap_db_dir` | `null` | Where ASTAP star tiles live. `null` → `CASSA_ASTAP_DB`, then `~/.cache/cassa-photometry/astap`. |
| `astap_db_series` | `auto` | `auto` picks from the frame's own field height; or name `d50`, `d05`, `g05`. |
| `astap_db_url` | SourceForge `star_databases/` | Where the published archives live. Only ZIP-published series can be fetched incrementally. |
| `astap_db_download` | `true` | `false` selects but never fetches, and reports which tiles are missing. |
| `astap_tile_neighbours` | `1` | Tiles either side of the field centre. `1` covers a field on a tile edge. |
| `astrometry_index_dir` | `null` | A full local astrometry.net index set. `null` → `CASSA_ASTROMETRY_INDEX`, then `./astrometry_data`. A populated directory wins over the on-demand cache. |
| `index_url` | `https://data.astrometry.net/` | Base URL of the index series. |
| `index_token` | `null` | Bearer token, for a private mirror. |
| `index_cache_dir` | `null` | `null` → `CASSA_INDEX_CACHE`, then `~/.cache/cassa-photometry/astrometry`. |
| `index_cache_gb` | `20.0` | Cache size budget. |
| `index_search_radius_deg` | `3.0` | Search cone when selecting index tiles. |
| `index_blind_radius_deg` | `15.0` | The wider cone used by the fallback pass after a first solve fails. |
| `index_scale_lo_frac` | `0.30` | Quad-scale window floor, as a fraction of the field size. The finest indexes are by far the largest files and rarely what solves a frame, so the first pass asks for 30–100%. |
| `index_scale_lo_frac_wide` | `0.10` | The widened fallback pass drops to the full 10%. |
| `index_scale_hi_frac` | `1.00` | Quad-scale window ceiling. |
| `index_download` | `true` | `false` still selects, and reports the `cassa-index-fetch` command that would get each missing file. |
| `solve_scale_tolerance` | `0.25` | Half-width of the pixel-scale window handed to the solver, as a fraction. Deliberately generous: a wrong stated scale makes a solve **fail**, where no scale only makes it slow. |
| `solve_scale_tolerance_wide` | `1.0` | The fallback pass's window. |
| `solve_scale_warn_frac` | `0.05` | Warn when the solved plate scale disagrees with the header's claim by more than this fraction. This is how a wrong `SECPIX` gets noticed. |
| `solve_timeout_s` | `120.0` | First-pass solve timeout. |
| `solve_timeout_wide_s` | `600.0` | Fallback-pass timeout. |

### `phase3:` — photometry, zero point, catalogs

| Key | Default | Meaning |
|---|---|---|
| `steps` | all on | See [section 4](#4-the-step-registry). |
| `fwhm` | `3.5` | **Fallback only.** Apertures are sized from the measured `FWHMPX`; a fixed aperture whatever the seeing is worth up to 0.3 mag of per-epoch zero-point error. |
| `detection_threshold` | `5.0` | Detection threshold in sigma. |
| `detect_npixels` | `5` | Minimum connected pixels above threshold. |
| `deblend_nlevels` | `32` | Deblending levels. |
| `deblend_contrast` | `0.001` | Deblending contrast ratio. |
| `aperture_r_factor` | `2.0` | Aperture radius, in units of the measured FWHM. |
| `annulus_in_factor` | `5.0` | Background annulus inner radius, in FWHM. Must sit outside the PSF wings — at 3–4× FWHM a Moffat still puts ~1.25% of the star's light in the annulus. |
| `annulus_out_factor` | `8.0` | Annulus outer radius, in FWHM. |
| `zp_sigma_clip` | `3.0` | Sigma clip applied to the per-star zero points. |
| `zp_match_tol_arcsec` | `2.0` | Fixed cross-match radius, used when the WCS carries no `ASTRMS`. |
| `match_radius_sigma` | `3.0` | Cross-match radius = this × `ASTRMS`… |
| `match_radius_min_arcsec` | `1.0` | …clamped to this floor. The floor matters: reference positions carry their own error and stars have moved since the catalog epoch. |
| `match_radius_max_arcsec` | `5.0` | …and this cap. |
| `zp_reject_flagged` | `true` | Reject calibrators with a saturated/bad/cosmic-ray pixel in the aperture. The brightest catalog stars saturate first and are exactly the ones inverse-variance weighting trusts most. |
| `psf_photometry` | `true` | PSF-fitted photometry. *(Deprecated duplicate of `steps.psf_photometry`.)* |
| `aperture_correction` | `true` | Measure an aperture correction from a curve of growth, making the zero point a **total-flux** zero point. *(Deprecated duplicate of `steps.aperture_correction`.)* |
| `apcor_total_factor` | `5.0` | Radius treated as "total" in the curve of growth, in FWHM. |
| `ellipticity_star_max` | `0.15` | Objects rounder than this are treated as stars by `cassa-verify`. |
| `ab_flux_zero_jy` | `3631.0` | Zero-flux of the AB system, in Jy. |
| `apply_ab_offset` | `true` | Apply the band's AB−Vega offset in the Jy conversion. Without it, B fluxes are ~8% wrong. |
| `keep_negative_flux` | `true` | Keep non-detections in the catalog with a limiting magnitude, rather than dropping them and biasing faint number counts. |
| `catalog_cache` | `true` | Persist every successful reference-catalog query. |
| `catalog_cache_dir` | `null` | `null` → `CASSA_CATALOG_CACHE`, then `~/.cache/cassa-photometry/catalogs`. |
| `catalog_cache_days` | `180.0` | Refresh a cached query older than this. `0` keeps it forever. |
| `offline` | `false` | Never touch the network; use the cache only, and fail clearly when a field is not in it. |

### `phase4:` — diagnostics

Phase 4 measures rather than decides, so these only steer what it looks at.

| Key | Default | Meaning |
|---|---|---|
| `steps` | all on | See [section 4](#4-the-step-registry). |
| `fwhm_guess` | `3.5` | Starting FWHM for star detection. Shared with the FWHM measurement phases 1–2 depend on. |
| `detection_threshold` | `5.0` | Detection threshold in sigma. |
| `max_stars` | `25` | Stars fitted when estimating the PSF. |
| `cutout` | `15` | Cutout size in pixels for each fitted star. |
| `near_saturation_frac` | `0.95` | A pixel within this fraction of the frame maximum counts as near-saturated in the stage 0/1 map. Deliberately independent of the DQ `SATURATED` bit, so the two can be compared. |
| `limiting_snr` | `5.0` | SNR defining the limiting magnitude. |
| `mag_binsize` | `0.5` | Magnitude bin width for number counts. |
| `hist_bins` | `30` | Histogram bins. |

### Plate-solver availability by platform

Any one backend is enough. This is why ASTAP is tried first — it is the only one
published for every platform the pipeline runs on.

| Platform | `astap` | `solve-field` | `astrometry-py` |
|---|:---:|:---:|:---:|
| Linux x86-64 | ✓ | ✓ | ✓ |
| Linux aarch64 | ✓ | — | — |
| macOS Intel | ✓ | ✓ | ✓ |
| macOS Apple Silicon | ✓ | — | ✓ |
| Windows x64 / x86 / ARM64 | ✓ | — | — |

conda-forge publishes `astrometry` for `linux-64` and `osx-64` only; PyPI's
`astrometry` publishes no aarch64 and no Windows wheel. Verified 2026-09-09
against both indexes.

---

## 3. Environment variables

| Variable | Effect |
|---|---|
| `CASSA_ASTROMETRY_INDEX` | A full local astrometry.net index set. Wins over the on-demand cache; nothing is downloaded. |
| `CASSA_INDEX_CACHE` | Where on-demand index files are cached. Default `~/.cache/cassa-photometry/astrometry`. |
| `CASSA_INDEX_TOKEN` | Bearer token for a private index mirror. |
| `CASSA_ASTAP_DB` | Where ASTAP star tiles live. Default `~/.cache/cassa-photometry/astap`. |
| `CASSA_CATALOG_CACHE` | Where reference-catalog queries are cached. Default `~/.cache/cassa-photometry/catalogs`. |
| `XDG_CACHE_HOME` | Relocates all three caches at once, in the usual way. |

Resolution order is always **config value → environment variable → default**, so
a config file wins over the environment.

---

## 4. The step registry

Every phase declares its steps with dependency tokens. A step may be excluded
freely; it may be reordered anywhere its `requires` are still satisfied by an
earlier `provides`. An order that breaks a real dependency is **refused when the
config is read**, naming the violation. `movable: no` marks a step nested inside
another, with no seam to move it to.

### Phase 1 — instrument signature removal

| Step | Requires | Provides | Movable |
|---|---|---|---|
| `linearity` | — | `linearised` | yes |
| `overscan` | `linearised` | `trimmed` | yes |
| `bad_pixel_mask` | — | `bpm` | yes |
| `bias` | `trimmed` | `bias_subtracted` | yes |
| `dark` | `bias_subtracted` | `dark_subtracted` | yes |
| `flat` | `dark_subtracted` | `flat_fielded` | yes |
| `cosmic_rays` | `trimmed`, `bpm` | `cr_cleaned` | yes |
| `measure_fwhm` | `flat_fielded` | `frame_fwhm` | yes |

Cosmic-ray rejection before or after flat fielding is a genuine choice and is
permitted; flat before bias is not.

### Phase 2 — integration

| Step | Requires | Provides | Movable |
|---|---|---|---|
| `subtract_background` | — | `background_removed` | yes |
| `align` | — | `registered` | yes |
| `stack` | `registered`, `background_removed` | `stacked` | yes |
| `solve_wcs` | `stacked` | `wcs` | yes |
| `measure_fwhm` | `stacked` | `master_fwhm` | **no** |
| `visual_qa` | `stacked`, `wcs` | `qa_pdf` | yes |

### Phase 3 — photometry

| Step | Requires | Provides | Movable |
|---|---|---|---|
| `aperture_correction` | — | `apcor` | **no** |
| `zero_point` | — | `zero_point` | yes |
| `flux_calibration` | `zero_point` | `fluxcal` | yes |
| `catalog` | — | `catalog` | yes |
| `psf_photometry` | `catalog` | `psf_mag` | **no** |
| `classification` | `catalog`, `psf_mag` | `classes` | **no** |

### Phase 4 — diagnostics

| Step | Requires | Provides | Movable |
|---|---|---|---|
| `stage_0_raw` | — | `stage_0` | yes |
| `stage_1_calibrated` | — | `stage_1` | yes |
| `stage_2_master` | — | `stage_2` | yes |
| `stage_3_photometry` | — | `stage_3` | yes |

The four stages are independent, so any subset in any order is valid.

### Configuring a plan

```yaml
phase1:
  steps:
    cosmic_rays: false                    # a boolean, or...
    exclude: [cosmic_rays, flat]          # ...a list, which reads better in bulk
    order: [linearity, overscan, bad_pixel_mask, bias, dark,
            cosmic_rays, flat, measure_fwhm]
    custom: {my_step: local/my_steps.py:my_step}
```

`--show-plan` prints the resolved order for any phase without reducing
anything. See [`CUSTOMIZING.md`](CUSTOMIZING.md) for the full treatment.

---

## 5. FITS keywords

### Read from raw frames

The header is always the first authority; an instrument profile supplies only
what the header cannot say.

| Keyword | Used for |
|---|---|
| `IMAGETYP` | Frame classification (bias / dark / flat / science). Folder names are documentation, not classification. |
| `FILTER` | Filter identity, mapped to a stacking label and a science band. |
| `EXPTIME`, `EXPOSURE` | Exposure time, for dark matching and flux normalisation. |
| `EGAIN`, `GAIN` | System gain in e⁻/ADU. On CMOS, `GAIN` is often the unitless *setting* — write `EGAIN`. |
| `READNOIS` | Read noise in e⁻. |
| `READOUTM` | Readout mode (`HCG`/`LCG`). Without it, the mode is inferred from the gain setting. |
| `SATURATE`, `SATLEVEL`, `FULLWELL` | Saturation level in raw ADU. |
| `XBINNING`, `YBINNING` | Binning. A binning mismatch against the masters is a hard skip. |
| `XPIXSZ`, `YPIXSZ`, `FOCALLEN` | Plate scale, when `SECPIX`/`PIXSCALE` are absent. |
| `SECPIX`, `PIXSCALE` | Stated plate scale, used as a solver hint and checked against the solution. |
| `OBJCTRA`, `OBJCTDEC`, `TARGRA`, `TARGDEC` | Pointing, for the solver hint and the sky-data selection. |
| `DATE-OBS`, `MJD-OBS`, `TIMESYS` | Epoch binning. |
| `SITELAT`, `SITELONG` | Local-noon epoch binning and airmass. |
| `CCD-TEMP`, `SET-TEMP` | Calibration temperature matching. |
| `OBJECT`, `TARGNAME`, `OBJTYPE` | What is being observed, and therefore how it is measured. |
| `INSTRUME`, `TELESCOP` | Recorded through the reduction. |
| `CALSTAT` | Non-empty means the frame reports prior calibration; it is skipped unless `allow_precalibrated`. |
| `AIRMASS`, `GUIDERMS` | Carried through and reported. |

### Written by Phase 1

| Keyword | Meaning |
|---|---|
| `CALVERS` | Calibration vintage (currently `2`). Bumped whenever a change makes new output numerically incomparable with old, so mixed vintages are detectable rather than silently averaged. |
| `CALPLAN` | The steps that ran, in the order they ran. |
| `CALSKIP` | The steps that did not run. A partially-reduced frame says so. |
| `BUNIT` | `electron` after gain correction. |
| `GAINVAL` | The gain actually used, whatever supplied it. |
| `BIASSUB` | Whether the master dark had its bias pedestal removed. |
| `CRVETO` | Set when a cosmic-ray mask was discarded as implausible (`cr_max_fraction`). |
| `FWHMPX`, `FWHMASEC` | Measured FWHM in pixels and arcsec. Phase 3 sizes apertures from this. |
| `FWHMNSTR`, `FWHMELL` | Stars fitted, and the fitted ellipticity. |
| `RAWDIR`, `RAWFILE` | Where the frame came from, so the raw tree is findable later. |
| `NCOMBINE`, `NREJECT` | Frames combined into a master, and pixels rejected. |
| `EXPNORM` | Exposure normalisation applied. |

### Written by Phase 2

| Keyword | Meaning |
|---|---|
| `STACKCNT` | Frames in the stack. |
| `TOT_EXP` | Total integrated exposure. |
| `EXPMEAN` | Mean per-frame exposure. |
| `EPOCHKEY` | The epoch bin this master belongs to. |
| `BKGLEVEL`, `BKGND` | The background level removed before stacking. |
| `PIXSCALE` | Plate scale of the stack. |
| `STEPPLAN`, `STEPSKIP` | The resolved step plan, and what was skipped. |
| `WCSFROM` | Where the WCS came from. |
| `WCSSOLVR` | Which backend solved it — `astap`, `solve-field` or `astrometry-py`. |
| `ASTRMS` | Total astrometric RMS residual, in arcsec. |
| `CRDER1`, `CRDER2` | FITS-standard random error per axis, in degrees. |
| `ASTNSTAR` | Stars matched in the residual measurement. |
| `ASTRMSRC` | Which route produced `ASTRMS` — the solver's own matched-star table, or the reference catalog. This is what makes the backends equivalent. |
| `CTYPE1/2`, `CUNIT1/2`, `CRPIX1/2`, `CRVAL1/2`, `CD*_*` | The WCS itself. |

### Written by Phase 3

| Keyword | Meaning |
|---|---|
| `MAGZERO` | Filter-wise zero point. |
| `MAGZERR` | Its uncertainty. |
| `NZPSTARS` | Calibrator stars that survived clipping and flag rejection. |
| `PHOTREF` | Which reference catalog supplied them. |
| `PHOTSYS` | `Vega` or `AB` for this band. |
| `PHOTWAVE` | Effective wavelength, in nm. |
| `ABOFFSET` | The AB−Vega offset applied. |
| `ZPAPER` | The aperture the zero point was measured in. |
| `ZPAB` | The zero point on the AB system. |
| `ZPSCOPE` | Whether the zero point is an aperture or a total-flux one. |
| `APCOR` | Aperture correction, in magnitudes. |
| `APCORRMS` | Its scatter. |
| `APCORMTH` | How it was measured. |
| `FWHMUSED` | The FWHM apertures were actually sized from. |
| `FLUXCAL` | Set on `*_fluxcal.fits`; `BUNIT` becomes `Jy/pixel`. |
| `STEPPLAN`, `STEPSKIP` | The resolved step plan, and what was skipped. |

---

## 6. Catalog columns

`*_catalog.csv`, one row per detected source.

| Column | Meaning |
|---|---|
| `NUMBER` | Segmentation label. |
| `ALPHA_J2000`, `DELTA_J2000` | Sky position, degrees. |
| `X_IMAGE`, `Y_IMAGE` | Pixel position, **1-indexed** (FITS convention). |
| `FLUX_ISO`, `FLUXERR_ISO` | Isophotal flux and its error. |
| `MAG_INST`, `MAGERR_INST` | Instrumental magnitude, before the zero point. |
| `MAG_ISO`, `MAGERR_ISO` | Isophotal magnitude. **Not a total magnitude** — see below. |
| `MAGERR_STAT` | The statistical part of the error alone. |
| `MAG_APER`, `MAGERR_APER`, `FLUX_APER`, `FLUXERR_APER` | Measured in the zero point's own aperture. |
| `MAG_AUTO`, `MAGERR_AUTO`, `FLUX_AUTO`, `FLUXERR_AUTO` | Kron elliptical — the standard total magnitude. |
| `KRON_RADIUS` | The Kron radius used. |
| `MAG_PSF`, `MAGERR_PSF` | PSF-fitted — the matched filter for a point source. |
| `CHI2_PSF` | Goodness of the PSF fit. |
| `PSF_MINUS_AUTO` | What classification is based on. |
| `MAG_BEST`, `MAGERR_BEST` | **The default choice.** PSF for point sources, Kron for extended. |
| `CLASS` | `STAR`, `EXTENDED`, `AMBIGUOUS`, `SATURATED` or `EDGE`. Omitted when classification is off. |
| `CLASS_STAR` | Stellarity, 0–1. Omitted when classification is off. |
| `LIMIT_MAG` | Limiting magnitude at this position, for non-detections. |
| `SNR` | Signal-to-noise from the ERR plane. |
| `FLAGS` | DQ flags found inside the source's segment — see section 7. |
| `ISOAREA_IMAGE` | Segment area in pixels. |
| `FLUX_RADIUS` | Half-light radius. |
| `FWHM_IMAGE` | Per-source FWHM. |

`MAG_ISO` captures a brightness-dependent fraction of a source, so an
aperture-derived zero point applied to it produces a **tilt** rather than an
offset — against simulated truth it runs from +0.18 mag at V=12 to +1.70 mag at
V=16.5. It is kept under its own name for compatibility only.

When classification is off, `CLASS` and `CLASS_STAR` are **omitted** rather than
filled with a placeholder, so a downstream cut on either fails loudly instead of
silently selecting nothing. `MAG_BEST` falls back to Kron for every source.

---

## 7. Data quality (DQ) bits

The `DQ` extension is an integer bitmask; a pixel may carry several.

| Bit | Value | Name | Meaning |
|---|---|---|---|
| 0 | `1` | `SATURATED` | At or above the saturation level. |
| 1 | `2` | `BAD_PIXEL` | Hot or dead, from the bad-pixel mask. |
| 2 | `4` | `COSMIC_RAY` | Flagged by cosmic-ray rejection. |
| 3 | `8` | `NO_DATA` | NaN, or outside the warp footprint. |
| 4 | `16` | `REJECTED` | Some contributing frames were rejected here (stacks only). |

`0` means good. `cassa_photometry.fits_utils.build_dq()` composes them and
`DQ_FLAG_NAMES` maps them back to names.

---

## 8. Python API

Everything is importable; the CLI is a thin layer over these.

```python
from cassa_photometry.config import load_config
from cassa_photometry.instruments.registry import get_profile
from cassa_photometry.phase1_calibration import pipeline as phase1

config = load_config("my_config.yaml")
profile = get_profile(config.instrument, config)
phase1.run("raw/", "work/phase1", config=config, instrument=profile)
```

### Shared infrastructure

| Module | What it provides |
|---|---|
| `config` | `load_config(path)`, `PipelineConfig` and the per-phase dataclasses, `ConfigError`. |
| `steps` | `registry_for(phase)`, `resolve_plan(phase, toggles)`, `resolved_names()`, `describe_plan()`, `load_custom_step()`, `Step`, `StepPlanError`. |
| `fits_utils` | `write_mef()`, `read_mef()`, `build_dq()`, the `DQ_*` constants, `CALVERS`. |
| `paths` | `phase_dir()`, `sibling_phase_dir()`, `find_phase_dir()`, `find_raw_frames()`, `raw_tree_summary()`. |
| `psf` | `estimate_fwhm()`, `detect_stars()`, `build_epsf()`, `background_rms()`, `image_stats()`. |
| `logging_utils` | `get_logger()`, `set_default_level()`. |
| `targets` | `Target`, `TargetRegistry`, `normalise_type()`. |
| `doctor` | `run_checks()`, `report()`, and each `check_*` individually. |

### Instruments

| Module | What it provides |
|---|---|
| `instruments.base` | `InstrumentProfile` — the interface. `get_gain()`, `get_read_noise()`, `get_saturation()`, `get_pixel_scale()`, `get_image_type()`, `standardize_filter()`, `science_band()`, `apply_linearity()`, `get_overscan_region()`, `get_subframe_origin()`, `already_calibrated()`, and the `*_from_header()` pair each of them consults first. |
| `instruments.registry` | `get_profile(name, config)`, `available_profiles()`. |
| `instruments.cassa` | `Cassa8InchProfile`, including `set_detector_curve()` for loading measured PTC results. |
| `instruments.overrides` | `with_detector_overrides(profile, detector)` — layers a config `detector:` block over any profile. |

### Phase 1 — calibration

| Module | What it provides |
|---|---|
| `phase1_calibration.pipeline` | `run()`, `build_master_bias()`, `build_master_dark()`, `build_master_flat()`. |
| `phase1_calibration.processor` | `UniversalProcessor.process_science_frame()`, `CalibratedFrame`, `subtract_scaled_dark()`, `crop_master_to_frame()`, `implausible_cosmic_rays()`. |
| `phase1_calibration.data_models` | `StandardCCD`, `load_standardized_ccds()`. |

### Phase 2 — integration and astrometry

| Module | What it provides |
|---|---|
| `phase2_integration.pipeline` | `run()`, `IntegrationPipeline`. |
| `phase2_integration.wcs` | `WCSSolver.solve()` — solving, the rescue pass, and the residual record. |
| `phase2_integration.solvers` | `get_solver()`, `available_backends()`, `BACKENDS`. |
| `phase2_integration.solvers.base` | `Solver`, `SolveHints` (`widened()`, `without_scale()`), `SolveResult` (`residuals()`), `solved_pixel_scale()`. |
| `phase2_integration.solvers.astap` | `AstapSolver`, `find_binary()`. |
| `phase2_integration.solvers.solvefield` | `SolveFieldSolver`. |
| `phase2_integration.solvers.inprocess` | `InProcessSolver`. |
| `phase2_integration.astrometry_qc` | `measure()`, `pair_by_position()`, `detect_sources()`, `reference_positions()` — the backend-independent residual. |
| `phase2_integration.math_utils` | `MathEngine`: `register_plane()`, `register_variance()`, `register_mask()`, `weighted_stack()`, `calc_scale()`, `extract_2d_background()`. |
| `phase2_integration.io` | `FITSHandler`: `extract_metadata()`, `load_planes()`, `save_master()`. |
| `phase2_integration.epochs` | `epoch_key()`. |
| `phase2_integration.hardware` | `HardwareManager`: `allocate_resources()` honours `cpu_fraction`/`max_cores`; `estimate_resources()` returns the envelope as a dict (cores used vs available, peak RAM vs total, runtime, `fits_in_budget`); `print_estimations()` logs it; `describe()` renders the notebook summary. |
| `phase2_integration.visuals` | `VisualQAGenerator`, `find_raw_frame()`. |
| `phase2_integration.models` | `TargetGroup`. |

### Phase 3 — photometry

| Module | What it provides |
|---|---|
| `phase3_photometry.pipeline` | `run()`, `detect_band()`. |
| `phase3_photometry.engine` | `UniversalPhotometryEngine`: `calculate_local_zero_point()`, `generate_full_catalog()`, `export_flux_calibrated_image()`, `resolve_fwhm()`, `search_radius()`; `combine_zeropoints()`. |
| `phase3_photometry.catalogs` | `fetch_reference_catalog()` and the individual `query_gaia_synthetic()`, `query_apass()`, `query_panstarrs()`, `query_sdss()`. |
| `phase3_photometry.apcor` | `measure()`, `from_curve()`, `ApertureCorrection.at()`. |
| `phase3_photometry.morphology` | `classify()`, `Classification`. |
| `phase3_photometry.photsys` | `system_of()`, `ab_offset()`, `effective_wavelength_nm()`, `describe()`. |
| `phase3_photometry.verify` | `run()`, `verify_calibration()`, `match_radius_arcsec()`, `cross_match_and_report()`. |

### Phase 4 — diagnostics

| Module | What it provides |
|---|---|
| `phase4_diagnostics.pipeline` | `run()` — returns the metrics dict it also writes as `metrics.json`. |
| `phase4_diagnostics.stages` | `diagnose_raw()`, `diagnose_calibrated()`, `diagnose_master()`, `diagnose_photometry()`. |
| `phase4_diagnostics.plots` | `Report`, `zscale_panel()`, `histogram_panel()`, `scatter_panel()`, `radial_profile_panel()`, `text_panel()`, `new_page()` and friends. |
| `phase4_diagnostics.psf` | The PSF estimation used by every stage. |

### Sky data

| Module | What it provides |
|---|---|
| `astap_db` | `AstapTileStore.ensure()`, `required_tiles()`, `band_of()`, `tile_of()`, `series_for_fov()`, `resolve_db_dir()`, `AstapDatabaseError`. |
| `astrometry_index` | `required_indexes()`, `build_store()`, `IndexManifest`, `IndexEntry`, `LocalIndexStore`, `HttpIndexStore`, `field_size_arcmin()`. |
| `index_fetch` | `run()` — behind `cassa-index-fetch`. |

### Simulation

| Module | What it provides |
|---|---|
| `simulate` | `run(preset, outdir, seed)`. |
| `simulate.scene` | `Scene`: `render()`, `truth_table()`, `flux_e_per_s()`, `sky_e_per_s()`, `local_fwhm()`; `moffat_alpha()`, `enclosed_fraction()`. |
| `simulate.detector` | `DetectorModel`: `expose_science()`, `expose_bias()`, `expose_dark()`, `expose_flat()`. |
| `simulate.presets` | `NightModel.step()`, `delivered_fwhm_arcsec()`. |
| `simulate.frames` | `write_science()`, `write_calibration()`, `airmass_from_altitude()`. |

---

## 9. What `cassa-doctor` checks

One line per check, with the fix for anything that failed. The exit status is
the number of failures.

| Check | What it establishes |
|---|---|
| `platform` | Linux, macOS and Windows all report OK. Windows was a hard failure while no plate solver was published for it; ASTAP ended that and a native install has since been run. |
| `python` | Interpreter version against the floor in `pyproject.toml`. |
| `cassa-photometry` | The pipeline's own version, whether it is an editable install, and where it lives. |
| `dependencies` | Every declared runtime requirement, against its minimum. |
| `solver: astap` | Whether the ASTAP binary is on `PATH` or at `phase2.astap_path`. |
| `solver: solve-field` | Whether the Astrometry.net binary is available, and its version. On Windows it reports *not published for Windows* rather than *not on PATH* — there is nothing to go and install. |
| `solver: in-process` | Whether the importable `astrometry` is the PyPI solver rather than the conda bindings — they share a name, and only one of them can solve. |
| `solver: any` | Appears **only when none is usable**, and fails: phase 2 cannot solve a WCS, so phase 3 has no zero point. |
| `solver: in use` | Which backend a run would actually pick — usually what a support question is really about. |
| `astrometry indexes` | Whether a full local set is configured, or the on-demand cache will be used. |
| `index manifest` | Whether the shipped index catalogue is readable. |
| `index selection` | Whether the native HEALPix bindings are present, giving exact tile selection rather than the conservative fallback. |
| `index cache` | How much is cached, and where. |
| `astap database` | How many ASTAP star tiles are cached. Empty is fine — tiles are fetched per field. |
| `network: astrometry index server` | Reachability of `data.astrometry.net`. A **warning**, not a failure: the pipeline is expected to work offline from its caches. |
| `network: VizieR (APASS)` | Reachability of the reference-catalog service. Also a warning. |
| `work directory` | That the working directory can be written to. |

Skip the two network checks with `--no-network`.

Paste the whole output when asking for help — it names the versions, the
backend and the cache state in one place.
