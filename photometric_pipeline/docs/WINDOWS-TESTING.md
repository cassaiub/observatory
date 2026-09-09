# Verifying the native Windows install

Native Windows support is **provisional**: every piece is in place and
`install.ps1` builds a real environment, but no one has yet run it end to end on
Windows. This is the checklist that closes that gap. It takes about twenty
minutes.

Until it is done, `cassa-doctor` reports the platform as a warning rather than
OK, and the documentation recommends WSL for anything with a deadline.

## Why this is worth doing

Windows was unsupported because no plate solver was published for it. That
stopped being true:

| | Status |
|---|---|
| Plate solver | ASTAP publishes `win64`, `win32` and `win11_aarch64` command-line builds — each a single `astap_cli.exe` |
| Runtime dependencies | All 16 have a Windows wheel or are pure Python (checked against PyPI, 2026-09-09) |
| Pipeline code | No POSIX-only imports, no `fork`, no hardcoded path separators |

`solve-field` and the in-process solver remain unavailable on Windows —
conda-forge builds `astrometry` for `linux-64`/`osx-64` only, and PyPI's
`astrometry` ships no Windows wheel. **ASTAP is the only backend Windows has**,
which is the whole reason the pipeline gained it.

## What to run

On a Windows machine, in PowerShell (no admin needed):

```powershell
git clone https://github.com/cassaiub/observatory.git
cd observatory\photometric_pipeline
.\install.ps1
```

Then, in the same window:

```powershell
conda activate cassa-photometry     # or: .\.venv\Scripts\Activate.ps1
cassa-doctor
cassa --help
pytest -q
```

If you have a night of frames handy, the real test is a reduction:

```powershell
cassa-simulate --preset tiny --out sim
cassa-run -i sim\raw -o sim\work
```

## What to send back

The output of `cassa-doctor` is the single most useful thing — it names the
platform, the versions, which backend is in use, and the cache state in one
screen. Beyond that:

1. Whether `install.ps1` completed, and what it printed if it did not.
2. Whether `solver: astap` shows a path, and whether `solver: in use` says
   `astap`.
3. The `pytest -q` summary line.
4. If you ran a reduction: whether phase 2 wrote a WCS (`WCSSOLVR` and `ASTRMS`
   in a `Master_*.fits` header) and whether phase 3 produced a zero point.

## Things most likely to break, and what they look like

| Symptom | Cause |
|---|---|
| `astap.zip` is ~114 KB and `Expand-Archive` fails | SourceForge served the HTML interstitial. The installer sets a non-browser User-Agent to avoid this and checks the ZIP magic bytes, so it should refuse rather than proceed — if it *did* proceed, that is a bug worth reporting. |
| `solver: astap  not on PATH` after a successful install | The binary landed in `%LOCALAPPDATA%\cassa-photometry\bin` because the environment directory was not writable. Set `phase2.astap_path` to it, or add it to PATH. |
| `conda env create` fails | Usually a proxy or a long-path limit. `.\install.ps1 -Venv` avoids conda entirely; every dependency has a wheel. |
| A path-related error inside a phase | This is the interesting one — it would be a genuine portability bug, not an environment problem. Please send the traceback. |

## When it passes

Three things change, and they should change together:

1. `doctor.check_platform()` — Windows moves from `WARN` to `OK`, and the
   "provisional" wording goes.
2. `tests/test_platform_support.py` — `test_native_windows_is_provisional_not_a_failure`
   becomes an `OK` assertion.
3. The docs — the README, the manual's Part II, the participant handbook and
   `install.ps1`'s own header stop calling it provisional, and WSL becomes one
   supported route among several rather than the recommended one.

Delete this file at that point, or replace it with a line recording the date and
the Windows version it was verified on.
