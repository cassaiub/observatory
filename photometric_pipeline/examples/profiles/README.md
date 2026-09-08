# Instrument profiles

The pipeline ships profiles only for setups the CASSA Observatory operates
(`generic`, `cassa8`). Describing another setup never requires editing the
package.

| Your situation | What to do |
|---|---|
| The camera writes `EGAIN`, `READNOIS`, `SATURATE` | Nothing. Use `instrument: generic`. |
| It does not, but the values are constants | A `detector:` block in your config file. No code. |
| The values depend on the header (CMOS gain curve, multi-amp, odd filter names) | Copy `template.py`, then `instrument_module: my_profile.py:MyProfile`. |
| You want to distribute a profile | Register it under the `cassa_photometry.instruments` entry-point group. |

Files here:

- **`template.py`** — a commented starting point, with the rules a profile must
  respect and every optional override shown.
- **`itelescope.py`** — the retired iTelescope network profile, kept as a
  realistic worked example of a header-dependent profile. Read its docstring:
  two of its choices are there to be criticised rather than copied.

Run either without installing anything:

```bash
cassa-calibrate -i raw -o work/phase1 -c my_config.yaml
```

```yaml
# my_config.yaml
instrument_module: examples/profiles/itelescope.py:ITelescopeNetworkProfile
```

Full explanation, including how to keep a modified clone rebasable against
upstream, is in [`docs/CUSTOMIZING.md`](../../docs/CUSTOMIZING.md).
