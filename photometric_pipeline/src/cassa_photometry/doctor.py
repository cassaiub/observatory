"""Environment diagnosis: the first thing to run when something does not work.

``cassa-doctor`` answers, in one screen, the questions that otherwise cost a
round-trip with an off-campus user: is the package installed and from where, are
the scientific libraries new enough, is there a working plate solver, are the
astrometry indexes reachable, can the reference catalogs be queried, and is the
working directory writable.

Every check reports ``OK``, ``WARN`` or ``FAIL`` with a one-line reason, and a
``FAIL`` prints the command that fixes it. The exit status is non-zero only for
``FAIL``, so the script is usable in CI and at the end of ``install.sh``.
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import tempfile

OK, WARN, FAIL = "OK", "WARN", "FAIL"

_SYMBOL = {OK: "+", WARN: "!", FAIL: "x"}


class Check:
    """One diagnostic line: a status, what was checked, and what was found."""

    def __init__(self, name, status, detail, fix=None):
        self.name = name
        self.status = status
        self.detail = detail
        self.fix = fix

    def __repr__(self):  # pragma: no cover - debugging aid
        return f"Check({self.name!r}, {self.status!r}, {self.detail!r})"


def _version_of(module_name):
    """Installed version of a distribution, or None when it is absent."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(module_name)
    except PackageNotFoundError:
        return None


def check_platform():
    """Linux and macOS are verified; native Windows is provisional.

    Windows used to be reported as a hard failure, on the grounds that no plate
    solver was published for it. That is no longer true: ASTAP ships
    command-line builds for win64, win32 and ARM64, and every runtime
    dependency has a Windows wheel. So the pieces are all present, and
    ``install.ps1`` now installs them.

    What is still missing is evidence. No native Windows install has been
    verified end to end, so this reports ``WARN`` rather than ``OK`` -- enough
    to say "this may work and nobody has checked", without the exit status
    claiming a failure that may not exist. WSL remains the route the project
    tests, and it is real x86-64 Linux, so every instruction applies unchanged
    inside it.

    When a Windows install has been confirmed, this becomes ``OK`` and the
    wording goes with it.
    """
    system = platform.system()
    detail = f"{system} {platform.machine()}"
    if system == "Windows":
        return Check(
            "platform", WARN,
            f"{detail} (native Windows support is provisional -- not yet verified)",
            "It should work: run `.\\install.ps1`. If anything fails, "
            "`.\\install.ps1 -Wsl` prints the tested WSL route. Either way, "
            "please report what this command printed.",
        )
    return Check("platform", OK, detail)


def check_python():
    v = sys.version_info
    detail = f"{platform.python_version()} ({platform.system()} {platform.machine()})"
    if v < (3, 10):
        return Check("python", FAIL, detail, "cassa-photometry requires Python >= 3.10.")
    return Check("python", OK, detail)


def check_package():
    """The package itself: version, and whether it is an editable install."""
    from cassa_photometry import __version__

    try:
        import cassa_photometry

        location = os.path.dirname(os.path.dirname(os.path.abspath(cassa_photometry.__file__)))
    except Exception:  # pragma: no cover - import already succeeded to get here
        location = "unknown"
    editable = "editable" if location.endswith("src") else "installed"
    return Check("cassa-photometry", OK, f"{__version__} ({editable}) at {location}")


def check_dependencies():
    """Every declared runtime requirement, against the floor in pyproject.toml.

    The requirement list is read from the installed metadata rather than
    duplicated here, so this check cannot drift from the packaging.
    """
    from importlib.metadata import requires

    try:
        from packaging.requirements import Requirement
    except ImportError:  # pragma: no cover - packaging ships with pip
        return [Check("dependencies", WARN, "packaging is unavailable; versions unchecked")]

    try:
        declared = requires("cassa-photometry") or []
    except Exception:
        return [Check("dependencies", WARN, "package metadata not found; run pip install -e .")]

    checks = []
    for spec in declared:
        req = Requirement(spec)
        if req.marker is not None and not req.marker.evaluate():
            continue  # an extra, or a platform this does not apply to
        found = _version_of(req.name)
        if found is None:
            checks.append(
                Check(req.name, FAIL, "not installed", f'pip install "{spec}"')
            )
        elif req.specifier and found not in req.specifier:
            checks.append(
                Check(
                    req.name,
                    FAIL,
                    f"{found} does not satisfy {req.specifier}",
                    f'pip install -U "{spec}"',
                )
            )
        else:
            checks.append(Check(req.name, OK, found))
    return checks


def _solve_field_version(executable):
    try:
        out = subprocess.run(
            [executable, "--version"], capture_output=True, text=True, timeout=20
        )
        return (out.stdout or out.stderr).strip().splitlines()[0][:60]
    except Exception:
        return "version unknown"


def check_astap_database(config=None):
    """How many ASTAP star tiles are cached locally.

    Unlike the astrometry.net index set there is nothing to pre-download: tiles
    are fetched per field, a few megabytes at a time. So an empty cache is
    normal on a fresh install and is reported as information, not a problem.
    """
    from cassa_photometry.astap_db import resolve_db_dir

    if config is None:
        from cassa_photometry.config import load_config

        config = load_config()

    directory = resolve_db_dir(config)
    if not os.path.isdir(directory):
        return Check("astap database", OK,
                     f"{directory} (empty; tiles are fetched per field)")
    tiles = [n for n in os.listdir(directory) if not n.endswith(".part")]
    size = sum(os.path.getsize(os.path.join(directory, n)) for n in tiles) / 1048576.0
    if not tiles:
        return Check("astap database", OK,
                     f"{directory} (empty; tiles are fetched per field)")
    return Check("astap database", OK,
                 f"{len(tiles)} tile(s), {size:.0f} MB in {directory}")


def check_solvers():
    """Which plate-solving backends are usable.

    Three are supported and any one is sufficient. Note that conda-forge's
    ``astrometry`` package (which supplies ``solve-field``) installs a Python
    module *also* called ``astrometry``, distinct from the PyPI package of that
    name that provides the in-process solver -- so presence of the module says
    nothing until we look for ``Solver`` on it.
    """
    checks = []

    # ASTAP first: it is what `auto` prefers, and the one most likely to be
    # present on a machine where neither of the others can be installed.
    from cassa_photometry.phase2_integration.solvers.astap import find_binary

    astap = find_binary()
    if astap:
        checks.append(Check("solver: astap", OK, astap))
    else:
        checks.append(Check(
            "solver: astap", WARN, "not on PATH",
            ".\\install.ps1   # or unpack https://www.hnsky.org/astap.htm yourself"
            if platform.system() == "Windows" else
            "apt install astap-cli   # or https://www.hnsky.org/astap.htm",
        ))

    executable = shutil.which("solve-field")
    if executable:
        checks.append(
            Check("solver: solve-field", OK, f"{executable} ({_solve_field_version(executable)})")
        )
    else:
        checks.append(Check(
            "solver: solve-field", WARN,
            "not published for Windows" if platform.system() == "Windows"
            else "not on PATH"))

    try:
        import astrometry as _astrometry

        if hasattr(_astrometry, "Solver"):
            checks.append(
                Check(
                    "solver: in-process",
                    OK,
                    f"astrometry {_version_of('astrometry') or '?'} (PyPI solver)",
                )
            )
        else:
            checks.append(
                Check(
                    "solver: in-process",
                    WARN,
                    "the importable 'astrometry' is the astrometry.net binding, "
                    "not the PyPI solver",
                )
            )
    except ImportError:
        checks.append(Check(
            "solver: in-process", WARN,
            "PyPI 'astrometry' publishes no Windows wheel"
            if platform.system() == "Windows"
            else "PyPI 'astrometry' not installed"))

    if all(c.status != OK for c in checks):
        checks.append(
            Check(
                "solver: any",
                FAIL,
                "no usable plate solver -- phase 2 cannot solve a WCS",
                ".\\install.ps1   # installs ASTAP, the one backend Windows has"
                if platform.system() == "Windows" else
                "apt install astap-cli   "
                '# or: pip install "cassa-photometry[solver]"',
            )
        )
    else:
        # Which one a run would actually use, which is what a support question
        # is usually really about.
        try:
            from cassa_photometry.config import load_config
            from cassa_photometry.logging_utils import get_logger
            from cassa_photometry.phase2_integration.solvers import get_solver

            active = get_solver(get_logger("cassa_doctor"), load_config())
            if active is not None:
                checks.append(Check("solver: in use", OK, active.name))
        except Exception as exc:  # pragma: no cover - defensive
            checks.append(Check("solver: in use", WARN, f"could not resolve ({exc})"))
    return checks


def check_astrometry_indexes(config=None):
    """Where index files will come from: a local directory, or the cache."""
    from cassa_photometry.config import load_config

    config = config or load_config()
    checks = []

    index_dir = config.resolve_astrometry_index_dir()
    if index_dir and os.path.isdir(index_dir):
        files = [f for f in os.listdir(index_dir) if f.startswith("index-") and f.endswith(".fits")]
        total_gb = sum(
            os.path.getsize(os.path.join(index_dir, f)) for f in files
        ) / 1e9
        if files:
            checks.append(
                Check(
                    "astrometry indexes",
                    OK,
                    f"{len(files)} local file(s), {total_gb:.1f} GB in {index_dir}",
                )
            )
        else:
            checks.append(
                Check("astrometry indexes", WARN, f"{index_dir} exists but holds no index-*.fits")
            )
    else:
        checks.append(
            Check(
                "astrometry indexes",
                OK,
                f"no local set ({index_dir}); files are fetched on demand",
            )
        )

    # The manifest is what makes on-demand fetching possible at all.
    try:
        from cassa_photometry.astrometry_index import IndexManifest, healpix_available

        manifest = IndexManifest.default()
        if len(manifest):
            available_gb = sum(e.bytes for e in manifest.entries) / 1e9
            checks.append(
                Check("index manifest", OK,
                      f"{len(manifest)} files catalogued ({available_gb:.0f} GB available "
                      f"from {config.phase2.index_url})")
            )
        else:
            checks.append(
                Check("index manifest", WARN,
                      "not built; selection falls back to the local directory",
                      "python tools/build_index_manifest.py --from-server")
            )
        checks.append(
            Check("index selection", OK, "exact (HEALPix bindings present)")
            if healpix_available() else
            Check("index selection", WARN,
                  "geometric fallback (no HEALPix bindings) -- downloads a superset")
        )
    except Exception as exc:  # pragma: no cover - defensive
        checks.append(Check("index manifest", WARN, f"unreadable ({exc})"))

    cache = config.resolve_index_cache_dir()
    if os.path.isdir(cache):
        cached = [f for f in os.listdir(cache) if f.endswith(".fits")]
        size_gb = sum(os.path.getsize(os.path.join(cache, f)) for f in cached) / 1e9
        checks.append(Check("index cache", OK,
                            f"{len(cached)} file(s), {size_gb:.1f} GB in {cache}"))
    else:
        checks.append(Check("index cache", OK, f"empty ({cache})"))
    return checks


def check_network(timeout=8):
    """Reachability of the services the pipeline fetches from.

    A failure here is a WARN, not a FAIL: the pipeline is expected to work
    offline from its caches.
    """
    import urllib.request

    targets = [
        ("astrometry index server", "https://data.astrometry.net/"),
        ("VizieR (APASS)", "https://vizier.cds.unistra.fr/"),
    ]
    checks = []
    for name, url in targets:
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                checks.append(Check(f"network: {name}", OK, f"reachable (HTTP {response.status})"))
        except Exception as exc:
            checks.append(
                Check(
                    f"network: {name}",
                    WARN,
                    f"unreachable ({type(exc).__name__}); offline runs use the cache",
                )
            )
    return checks


def check_writable(path=None):
    path = path or os.getcwd()
    try:
        with tempfile.NamedTemporaryFile(dir=path):
            pass
        return Check("work directory", OK, f"{path} is writable")
    except Exception as exc:
        return Check(
            "work directory", FAIL, f"{path} is not writable ({type(exc).__name__})",
            "Run from a directory you own, or pass -o /path/you/own.",
        )


def run_checks(config=None, skip_network=False):
    """Every check, in report order."""
    checks = [check_platform(), check_python(), check_package()]
    checks += check_dependencies()
    checks += check_solvers()
    checks += check_astrometry_indexes(config)
    checks.append(check_astap_database(config))
    if not skip_network:
        checks += check_network()
    checks.append(check_writable())
    return checks


def report(checks, stream=None):
    """Print the checks as a table and return the number of failures."""
    stream = stream or sys.stdout
    width = max(len(c.name) for c in checks)
    print("cassa-doctor", file=stream)
    print("=" * (width + 40), file=stream)
    for check in checks:
        print(f"[{_SYMBOL[check.status]}] {check.name:<{width}}  {check.detail}", file=stream)

    failures = [c for c in checks if c.status == FAIL]
    warnings = [c for c in checks if c.status == WARN]
    print("=" * (width + 40), file=stream)
    print(
        f"{len(checks) - len(failures) - len(warnings)} ok, "
        f"{len(warnings)} warning(s), {len(failures)} failure(s)",
        file=stream,
    )
    if failures:
        print("\nTo fix:", file=stream)
        for check in failures:
            if check.fix:
                print(f"  {check.name}: {check.fix}", file=stream)
    return len(failures)


def main(argv=None):
    """Allow ``python -m cassa_photometry.doctor`` as well as ``cassa-doctor``.

    install.sh and CI use the module form because they know the interpreter they
    just installed into but not necessarily where that environment's scripts
    landed.
    """
    import argparse

    from cassa_photometry.config import load_config

    parser = argparse.ArgumentParser(prog="cassa-doctor")
    parser.add_argument("--no-network", action="store_true")
    parser.add_argument("-c", "--config", default=None)
    args = parser.parse_args(argv)

    return report(run_checks(config=load_config(args.config),
                             skip_network=args.no_network))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(1 if main() else 0)
