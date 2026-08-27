# cassa-camchar

CMOS **sensor characterization** for the CASSA Observatory (QHYCCD miniCAM8M /
Sony IMX585). From a directory of calibration frames it measures — via the
**Photon Transfer Curve** method — the system gain, read noise, full-well
capacity, dynamic range, dark current (vs temperature), dark linearity
(amp-glow), quantum efficiency, and filter transmission, and writes a results
table plus a 9-panel dashboard.

The capture procedures (lab setup, FITS headers, exposure/temperature sweeps)
are documented in **Part I** of the observatory master documentation; this
package is the analysis code (**Part V**).

## What was refined

- **Installable package** with two commands and a central YAML config (no more
  hardcoded paths/parameters).
- **Robust statistics:** all means/variances are computed on a central ROI with
  sigma-clipping (rejecting hot pixels and vignetted corners) instead of raw
  full-frame `np.mean`/`np.var`.
- **Persisted results:** `characterization_results.csv` + `.json` (previously
  only a PNG was produced).
- **Honest fallbacks:** the fabricated dark-linearity data is gone (insufficient
  data is flagged); simulated QE/filter curves are tagged `synthetic`.

## Installation

```bash
conda activate pipeline_env      # provides numpy/scipy/astropy/matplotlib
cd camera_characterization
pip install -e . --no-deps --no-build-isolation
```

Registers `cassa-camchar-analyze` and `cassa-camchar-sim`.

## Usage

```bash
# Characterize a data directory (bias/ flat/ dark/ spectra/ inside)
cassa-camchar-analyze data --outdir camchar_output

# Generate a synthetic campaign with a known sensor model (for testing/demos)
cassa-camchar-sim --out /tmp/camchar_sim
cassa-camchar-analyze /tmp/camchar_sim --outdir /tmp/camchar_out
```

Outputs (in the output directory):

| File | Contents |
|------|----------|
| `characterization_results.csv` | per-gain table: system gain, read noise, full well, dynamic range |
| `characterization_results.json` | full results incl. dark current, linearity, QE, filter metrics |
| `CASSA_Sensor_Characterization_Report.png` | the 9-panel dashboard |

## Data layout

```
data/
  bias/     bias_gain{G}_{n}.fit             (>= 2 per gain)
  flat/     flat_gain{G}_exp{t}s_{n}.fit     (2 per exposure per gain)
  dark/     dark_temp{T}c_exp{t}s.fit        (temperature + exposure sweeps)
  spectra/  qe_calibration.csv, filter_*.csv (optional; synthetic if absent)
```

Grouping is by FITS header (`GAIN`, `EXPTIME`, `CCD-TEMP`), not filename.

## Development

```bash
pip install -e ".[dev]"
pytest
```
