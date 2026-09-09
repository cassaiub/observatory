# Windows: what is verified, and what to check next

Native Windows is supported. This file records how that happened, what the
first real run found, and the one thing still worth confirming — so nobody has
to reconstruct it from commit messages.

## Status

| | |
|---|---|
| `install.ps1` on Windows + Miniconda | **Verified** — completes with no errors |
| Package imports, `cassa-doctor` runs | **Verified** |
| Phase 1 (calibration) | **Verified** — ran clean |
| Phase 2 (integration + WCS) | **Fixed, needs a re-run** — see below |
| Phases 3–4 | Not yet exercised on Windows |

## What the first run found

Two defects, both Windows-only, both invisible on Linux and macOS — which is
exactly why they survived until someone ran the pipeline on Windows.

### 1. Phase 2 could not write its own master back

```
PermissionError: [WinError 32] The process cannot access the file because it is
being used by another process:
'...\work\phase2\Master_NGC7331_B_Photo_..._5fr.fits'
```

Astropy memory-maps image data by default, and a memory-mapped array keeps the
file handle alive for as long as the array exists — even after the `HDUList` is
closed. Phase 2 is precisely the shape that trips this: read a master, solve its
WCS, write it back over itself. POSIX lets you replace an open file; Windows
does not.

Fixed by routing every read in the package through `fits_utils.open_fits()`,
which turns mapping off, and by making `read_mef` return real copies rather than
views. `tests/test_windows_file_handles.py` pins the property on every platform,
because the symptom only appears on the one that is least likely to be in CI.

### 2. Index files were downloaded for a backend that never reads them

```
-> Solving with 8 index file(s); SolveHints(...)
-> ASTAP database: 10 tile(s) already cached.
```

Both lines for one solve, and only the second mattered. Phase 2 selected and
fetched Astrometry.net index files on every solve regardless of backend — about
246 MB per field handed to ASTAP, which ignores the argument entirely, on top of
the ~6 MB of star tiles it actually uses.

Fixed with a `uses_index_files` flag on the solver classes; ASTAP sets it
`False`. It also stops the "WCS solving will fail; run cassa-index-fetch"
warning, which was wrong and expensive advice for an ASTAP user.
`tests/test_solver_index_use.py` covers it.

## What to check next

A phase 2 re-run on Windows, to confirm the `PermissionError` is gone:

```powershell
conda activate cassa-photometry
cassa-run -i raw -o work --to 2
```

Then, if that passes, phases 3 and 4:

```powershell
cassa-photometry work\phase2
cassa-diagnose work\phase2 --raw raw
```

Useful to send back: the `cassa-doctor` output, whether a `Master_*.fits`
carries `WCSSOLVR` and `ASTRMS`, and whether phase 3 produced a zero point
(`MAGZERO`).

## Windows-specific notes

| Thing | Detail |
|---|---|
| Plate solver | **ASTAP only.** conda-forge builds `astrometry` for `linux-64`/`osx-64`; PyPI's `astrometry` ships no Windows wheel. This is why the pipeline gained ASTAP. |
| ASTAP download | SourceForge serves a *browser* an HTML "your download will start shortly" page instead of the file, and `Invoke-WebRequest` identifies as a browser by default — 114 KB of HTML named `astap.zip`. The installer asks with a non-browser user agent and verifies the ZIP magic bytes before trusting the result. |
| Binary location | `astap_cli.exe` goes into the environment's directory, or `%LOCALAPPDATA%\cassa-photometry\bin` if that is not writable. If `cassa-doctor` cannot find it, set `phase2.astap_path`. |
| Long paths | `conda env create` can fail on the 260-character limit. `.\install.ps1 -Venv` avoids conda entirely — every dependency has a Windows wheel. |
| Where to work | Anywhere. If you use WSL instead, keep the repository in your WSL home; reducing across `/mnt/c` is several times slower. |
