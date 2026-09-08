# Customizing the pipeline

Most changes people want do not need a code change, and the ones that do should
not need you to maintain a permanent fork. This document goes through the four
levels, cheapest first, and ends with how to keep a modified clone mergeable
with upstream.

| You want to… | Level |
|---|---|
| Turn a reduction step off | **1. Config** |
| Run the steps in a different order | **1. Config** |
| Add a step of your own | **1. Config** |
| Change a threshold, radius or tolerance | **1. Config** |
| Describe a camera the pipeline ships no profile for | **2. Instrument** |
| Change how a value is derived from the header | **2. Instrument** |
| Add a measurement, or change an algorithm | **3. Code** |
| Distribute your changes to others | **4. Package** |

---

## 1. Config: turn steps off, reorder them, add your own

Every tunable value lives in `cassa_photometry.config` and can be overridden
from a YAML file passed with `-c/--config`. Nothing else is needed.

```bash
cassa-run -i raw -o work --instrument cassa8 -c my_config.yaml
```

### Skipping a step

Each phase has a `steps:` block. Setting one to `false` skips it.

```yaml
phase1:
  steps:
    cosmic_rays: false     # short exposures, or you clean them elsewhere
    flat: false            # frames were flat-fielded before they reached you
phase2:
  steps:
    solve_wcs: false       # stack only; solve later, or the field has no solution
phase3:
  steps:
    psf_photometry: false  # faster; MAG_BEST then falls back to Kron
```

The available steps are:

| Phase | Steps |
|---|---|
| `phase1` | `overscan`, `bias`, `dark`, `flat`, `linearity`, `cosmic_rays`, `bad_pixel_mask`, `measure_fwhm` |
| `phase2` | `align`, `subtract_background`, `stack`, `solve_wcs`, `measure_fwhm`, `visual_qa` |
| `phase3` | `zero_point`, `aperture_correction`, `flux_calibration`, `catalog`, `psf_photometry`, `classification` |
| `phase4` | `stage_0_raw`, `stage_1_calibrated`, `stage_2_master`, `stage_3_photometry` |

A few of these change *what is produced*, not merely whether a correction is
applied, so they are worth stating explicitly:

| Step | What "off" produces |
|---|---|
| `phase2.align` | Frames are combined on the pixel grid they already share (a tracked mount, no dither). A frame of a different shape is an error rather than being resampled. No star matching means no transparency normalisation, so the flux scale stays 1.0. |
| `phase2.subtract_background` | The sky is measured (`BKGLEVEL` is still written) but kept in the data. This is the setting for extended sources — a background model follows the outer disk of a resolved galaxy and eats it. |
| `phase2.stack` | Aligned per-frame products (`Aligned_<group>_NNNofM.fits`, SCI/ERR/DQ) instead of one master, so the time axis survives into phase 3. Only the anchor is plate-solved; the others share its WCS, which is valid because registration put them on its grid. `min_frames_to_stack` does not apply. |
| `phase2.measure_fwhm` | No `FWHMPX` on the master — including no inherited one from the anchor frame. Phase 3 falls back to `phase3.fwhm`, so apertures are no longer matched to the seeing. |
| `phase3.flux_calibration` | The zero point is still measured and recorded; the full-size `_fluxcal.fits` is not written. |
| `phase3.classification` | `CLASS` and `CLASS_STAR` are omitted rather than filled with a placeholder, and `MAG_BEST` falls back to Kron (`MAG_AUTO`) for every source. A cut on `CLASS` then fails loudly instead of silently selecting nothing. |

### Reordering steps

`order:` sets an explicit sequence. It must list exactly the steps that run,
once each — `--show-plan` prints the list to start from.

```yaml
phase1:
  steps:
    # Reject cosmic rays BEFORE flat fielding rather than after.
    order: [linearity, overscan, bad_pixel_mask, bias, dark,
            cosmic_rays, flat, measure_fwhm]
```

**An order that breaks a dependency is refused when the config is read**, not
three hours into a reduction:

```
error: phase1.steps.order: 'flat' is placed before 'dark', but it needs what
that step produces ('dark_subtracted'). Move 'dark' earlier, or exclude it if
you do not want it at all.
```

This is deliberate. Order is not a free choice: overscan must come first because
it changes the array shape, and bias → dark → flat is physics, not preference.
Flat-fielding before bias subtraction produces a plausible, *silently wrong*
image — the same class of failure `CALSKIP`, `CRVETO` and the pre-calibrated
guard exist to prevent.

Equally, the pipeline does not over-constrain. Where two orders are both
defensible, no dependency links them and either is accepted:

| Genuine choice | Why both are defensible |
|---|---|
| `cosmic_rays` before or after `flat` | Detection behaves differently on a flat-fielded frame; neither is obviously right |
| `subtract_background` before or after `align` | Fitting before avoids interpolating the sky; fitting after measures it on the grid the stack uses |
| `aperture_correction` before or after `zero_point` | The correction can be measured first or applied to a combined ZP |

**Excluding is always legal**, even for a step others depend on: if nothing
provides a token, nothing requires it. `--skip bias` leaves `dark` and `flat`
running, which is what reducing pre-reduced frames needs.

**Some steps cannot be reordered at all.** `phase3.psf_photometry`,
`phase3.classification` and `phase3.aperture_correction` run *inside* another
step rather than beside it, as does `phase2.measure_fwhm`. They can be switched
off, but an `order` that moves them is refused rather than silently ignored.

### Adding your own step

```yaml
phase3:
  steps:
    custom:
      write_bright_list: local/my_steps.py:write_bright_list
    order: [aperture_correction, zero_point, flux_calibration, catalog,
            psf_photometry, classification, write_bright_list]
```

The step is a function taking keyword arguments only, so declare what you use
and swallow the rest with `**_`:

```python
def write_bright_list(*, engine, outdir, base, logger, **_):
    ...

write_bright_list.requires = ("catalog",)   # joins the dependency graph
```

Phase 3 passes `engine`, `path`, `outdir`, `base`, `logger`; phase 1 passes
`ccd`, `meta`, `state`, `logger`. Setting `requires` / `provides` makes a bad
placement a load-time error instead of an empty output file. See
[`examples/steps/my_steps.py`](../examples/steps/my_steps.py) for worked
examples. Your file is one upstream does not have, so `git pull` keeps merging.

### From the command line

Every phase command takes `--skip`, `--only` and `--show-plan`; they override
the config file.

```bash
cassa-calibrate -i raw -o work --skip cosmic_rays --skip flat
cassa-integrate work/phase1 --only stack
cassa-run -i raw -o work --show-plan          # print the plan, reduce nothing
cassa-run -i raw -o work --from 2 --to 3      # resume; do not recalibrate
```

On `cassa-run`, a step name that exists in two phases (`measure_fwhm`) must be
qualified as `phase2.measure_fwhm` — an ambiguous name is an error, not a guess.

### What was run is recorded

Every product states its own plan, so a customised reduction can never be
mistaken for a default one afterwards:

| Product | Cards |
|---|---|
| Phase 1 calibrated frames | `CALPLAN` (steps run, in order), `CALSKIP` |
| Phase 2 masters | `STEPPLAN`, `STEPSKIP` |
| Phase 3 `_fluxcal.fits` | `STEPPLAN`, `STEPSKIP` |
| Phase 4 `metrics.json` | `step_plans` for all four phases |

### Deprecated duplicates

Four settings used to exist twice and were combined with AND, so switching one
off left the other still blocking. The `steps` toggle is now the single
authority. The old names still work and still disable their step, but warn:

`phase1.apply_linearity`, `phase2.subtract_background`,
`phase3.psf_photometry`, `phase3.aperture_correction`.

**A skipped step is recorded.** Phase 1 writes `CALSKIP` into every calibrated
frame naming what it did not do, and each phase logs it. That is deliberate: a
partially-reduced frame that looks finished is how a wrong result gets
published.

### Common overrides

```yaml
instrument: cassa8
log_level: DEBUG

phase1:
  saturation_adu: 60000          # fallback only; the profile wins if it knows
  master_sigma_clip: 3.0         # 0 disables rejection in calibration stacks
  calibration_temp_tolerance_c: 3.0

phase2:
  solver: auto                   # auto | solve-field | astrometry-py
  epoch_bin: night               # none | night | "6h"
  stack_weight: point_source     # point_source | extended
  fwhm_reject_factor: 1.6        # 0 keeps every frame however soft
  index_download: true           # false: select but never fetch

phase3:
  offline: false                 # true: reference catalogs from cache only
  aperture_r_factor: 2.0         # aperture radius, in units of the measured FWHM
  zp_reject_flagged: true        # drop saturated/flagged calibrators
  keep_negative_flux: true       # keep non-detections, with a limiting magnitude
```

Every key is validated: a wrong type or a mistyped name is reported against its
full dotted path rather than silently ignored.

### Where the caches live

```bash
~/.cache/cassa-photometry/astrometry   # index files, fetched on demand
~/.cache/cassa-photometry/catalogs     # reference-catalog queries
```

Override with `phase2.index_cache_dir`, `phase3.catalog_cache_dir`, or the
`CASSA_INDEX_CACHE` / `CASSA_CATALOG_CACHE` environment variables. Reducing a
field once populates both, after which `--offline` works.

---

## 2. Instrument: describe your setup

See [`examples/profiles/`](../examples/profiles/) for the full worked examples.
In short:

**A `detector:` block** covers the common case — a camera that does not record
its own constants:

```yaml
instrument: generic
detector:
  gain: 1.4                 # e-/ADU
  read_noise: 7.0           # e-
  saturation_adu: 60000
  pixel_scale_arcsec: 0.40  # unbinned
  filter_map: {Sloan-R: R_Photo}
  science_bands: {Sloan-R: R}
```

These fill in **where the header is silent**. The frame stays the first
authority on its own data; `override_header: true` inverts that and logs every
card it discards.

**A profile class** is for anything header-dependent — a CMOS conversion-gain
curve, a multi-amplifier readout, an unusual filter wheel:

```yaml
instrument_module: my_profile.py:MyObservatoryProfile
```

No installation and no edit to this repository, which is what keeps your clone
rebasable.

---

## 3. Code: change an algorithm

If you have to change code, change it at a seam. These are the places designed
to be replaced, and the ones least likely to move under you:

| Seam | What it decides | Where |
|---|---|---|
| `InstrumentProfile` | Everything detector- and telescope-specific | `instruments/base.py` |
| `Solver` | How a WCS is obtained | `phase2_integration/solvers/base.py` |
| `IndexManifest` / `IndexStore` | Which astrometry indexes are used and where they come from | `astrometry_index.py` |
| `fetch_reference_catalog` | Which reference catalog calibrates a band | `phase3_photometry/catalogs.py` |
| `photsys` | What system a band is on | `phase3_photometry/photsys.py` |
| `apcor.measure` | How the aperture correction is derived | `phase3_photometry/apcor.py` |
| `morphology.classify` | How sources are classified | `phase3_photometry/morphology.py` |
| `TargetRegistry` | What is being observed | `targets.py` |
| `UniversalProcessor._detect_cosmics` | How cosmic rays are found and repaired | `phase1_calibration/processor.py` |

### When cosmic-ray rejection removes stars instead

astroscrappy finds cosmic rays as sharp positive features that its
fine-structure model cannot explain. In a crowded field -- a globular cluster,
a large resolved galaxy -- a great many stellar peaks answer that description,
and the step does not merely flag them: it *replaces* them. The result is a
frame with real starlight deleted, DQ flags sitting on the surviving stars, and
no sign of either in the calibrated product.

The pipeline will not let that pass silently. After each frame it measures what
share of the pixels the mask claims, and discards the mask whole if that share
exceeds `phase1.cr_max_fraction` (1% by default):

```
[WARNING] Cosmic-ray rejection flagged 11.4% of a frame, far more than any
          cosmic-ray rate can produce. The mask is being discarded ...
```

The threshold is physical rather than tuned. Cosmic rays arrive at a rate set by
the sky, not by what the telescope points at -- a few events per cm^2 per
minute, each marking a handful of pixels -- so even a large sensor in a long
exposure stays well under 0.1% of its pixels. The vetoed frame records `CRVETO`
in its header with the fraction that was rejected, so a frame that skipped the
step can never be mistaken for one that passed it.

The mask is discarded rather than merely left unapplied, because keeping the
flags is not the cautious half of the choice: phase 3 disqualifies any zero-point
calibrator whose aperture contains a flagged pixel, so a mask that lies on the
stars removes the very sources the zero point is measured from.

When this fires, reject cosmic rays in the phase 2 stack instead -- they do not
repeat between dithered frames, so a median across five of them removes what a
single frame cannot distinguish. To accept the mask anyway, raise the limit:

```yaml
phase1:
  cr_max_fraction: 1.0     # accept any mask, however large
```

### A worked example: two ways to remove cosmic-ray rejection

**By config** — the answer if you simply do not want the step:

```yaml
phase1:
  steps:
    cosmic_rays: false
```

**By code** — the answer if you want a *different* algorithm. Subclass the
processor and override the one method, rather than editing it in place:

```python
# local/my_processor.py
from cassa_photometry.phase1_calibration.processor import UniversalProcessor


class MyProcessor(UniversalProcessor):
    """Rejects cosmic rays with my own detector, not astroscrappy."""

    def _detect_cosmics(self, data, meta, inmask, ccd):
        mask = my_detector(data)
        return mask, my_repair(data, mask)
```

The point of the second form is that `git pull` still merges cleanly: your code
is in a file upstream does not have.

### Adding a catalog column

`generate_full_catalog` in `phase3_photometry/engine.py` builds an
`astropy.table.Table` named `out`. Adding a column is one line there. Keep the
SExtractor naming convention (`MAG_*`, `FLUX_*`, `*_IMAGE`) — anyone receiving
the catalog will expect it.

---

## 4. Package: distribute your changes

To share an instrument profile without anyone editing this repository, register
it under the entry-point group in your own `pyproject.toml`:

```toml
[project.entry-points."cassa_photometry.instruments"]
my_scope = "my_package.profiles:MyScopeProfile"
```

It then appears in `cassa-calibrate --instrument` for anyone who installs your
package. A plugin that fails to import is reported and skipped, never fatal.

---

## Keeping a modified clone mergeable

The single most useful habit: **put your changes in files upstream does not
have**, and reach them from config.

```bash
git remote add upstream <repo-url>
git fetch upstream
git rebase upstream/main          # or: git merge upstream/main
```

```
your-clone/
  local/                 # your code — upstream never touches this
    my_profile.py
    my_processor.py
  configs/
    cassa8_2026.yaml     # your settings
  src/cassa_photometry/  # leave alone wherever you can
```

If you must edit a shipped file, keep the edit small and in one place, and add a
comment saying why — it is what makes the conflict resolvable a year later when
you have forgotten.

### Before you fork, check it is not already configurable

```bash
python -c "
from cassa_photometry.config import PipelineConfig
import dataclasses, pprint
pprint.pp(dataclasses.asdict(PipelineConfig()))
"
```

That prints every setting the pipeline has, with its default.

---

## Checking your changes did what you meant

The pipeline ships a simulator with known truth, which is the honest way to
check that a modification improved something rather than merely changed it:

```bash
cassa-simulate --preset tiny --out sim
cassa-run -i sim/raw -o sim/work --instrument cassa8 -c my_config.yaml
```

`sim/truth_sources.csv` holds every source's true magnitude, and
`sim/truth_frames.csv` every frame's true seeing, transparency and zero point.
Compare your catalog against them before and after your change.

And run the tests — they encode a good deal of what the pipeline is supposed to
guarantee:

```bash
pytest -m "not network"
```
