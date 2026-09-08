#!/usr/bin/env bash
#
# cassa-photometry: clone and run this. Nothing else.
#
#   git clone <repo-url>
#   cd photometric_pipeline
#   ./install.sh
#
# It picks an environment, installs every dependency including a working plate
# solver, installs the pipeline, and then runs `cassa-doctor` to prove it works.
# Running it twice is safe -- it updates in place.
#
# Options:
#   --conda            force the conda route (needed on Linux aarch64)
#   --venv             force a plain python3 -m venv in ./.venv
#   --new-env          create a dedicated env even if one is already active
#   --env-name NAME    conda environment name (default: cassa-photometry)
#   --no-dev           skip pytest/ruff
#   --no-doctor        skip the closing health check
#   -h, --help         show this
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

MODE=""                       # conda | venv | current
ENV_NAME="cassa-photometry"
EXTRAS="[dev]"
RUN_DOCTOR=1
FORCE_NEW_ENV=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --conda)      MODE="conda" ;;
        --venv)       MODE="venv" ;;
        --new-env)    FORCE_NEW_ENV=1 ;;
        --env-name)   ENV_NAME="$2"; shift ;;
        --no-dev)     EXTRAS="" ;;
        --no-doctor)  RUN_DOCTOR=0 ;;
        -h|--help)    awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' "$0"; exit 0 ;;
        *) echo "install.sh: unknown option '$1' (try --help)" >&2; exit 2 ;;
    esac
    shift
done

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[33m[!] %s\033[0m\n' "$*"; }
die()  { printf '\033[31m[x] %s\033[0m\n' "$*" >&2; exit 1; }

# --- choose how to install ---------------------------------------------------
#
# Preference order, and why:
#   1. an already-active conda env  -- the user has told us where they work
#   2. conda, if available          -- it can supply solve-field, pip cannot
#   3. a venv                       -- always works; uses the in-process solver
if [[ -z "$MODE" ]]; then
    if [[ $FORCE_NEW_ENV -eq 0 && -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
        MODE="current"
    elif command -v conda >/dev/null 2>&1; then
        MODE="conda"
    else
        MODE="venv"
    fi
fi

PY=""       # the interpreter to install into
RUNNER=()   # prefix for running commands in that environment

case "$MODE" in
current)
    say "Using the active conda environment: $CONDA_DEFAULT_ENV"
    echo "    (pass --new-env to build a dedicated one instead)"
    PY="$(command -v python)"
    ;;
conda)
    command -v conda >/dev/null 2>&1 || die "conda not found. Install Miniforge, or re-run with --venv."
    if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        say "Updating conda environment '$ENV_NAME'"
        conda env update -n "$ENV_NAME" -f environment.yml --prune
    else
        say "Creating conda environment '$ENV_NAME'"
        conda env create -n "$ENV_NAME" -f environment.yml
    fi
    RUNNER=(conda run --no-capture-output -n "$ENV_NAME")
    PY="$("${RUNNER[@]}" python -c 'import sys; print(sys.executable)')"
    ;;
venv)
    command -v python3 >/dev/null 2>&1 || die "python3 not found."
    if [[ ! -d .venv ]]; then
        say "Creating virtual environment in ./.venv"
        python3 -m venv .venv
    else
        say "Reusing the virtual environment in ./.venv"
    fi
    PY="$REPO_ROOT/.venv/bin/python"
    "$PY" -m pip install --quiet --upgrade pip setuptools wheel
    ;;
esac

[[ -n "$PY" && -x "$PY" ]] || die "could not locate the environment's python."

# --- install the package -----------------------------------------------------
say "Installing cassa-photometry (editable)"
"$PY" -m pip install -e ".${EXTRAS}"

# --- make sure there is a plate solver ---------------------------------------
#
# Phase 2 needs one of two backends. `solve-field` (the Astrometry.net binary)
# cannot be pip-installed, so when it is missing we add the PyPI `astrometry`
# package, which solves in-process.
#
# These two MUST NOT coexist: conda-forge's astrometry package installs Python
# bindings that import under the same name as the PyPI solver, and whichever
# lands second wins. So this only ever runs when solve-field is absent.
say "Checking for a plate solver"
ENV_BIN="$(dirname "$PY")"
UNAME_S="$(uname -s 2>/dev/null || echo unknown)"
UNAME_M="$(uname -m 2>/dev/null || echo unknown)"

has_solve_field() {
    command -v "$ENV_BIN/solve-field" >/dev/null 2>&1 || command -v solve-field >/dev/null 2>&1
}

if has_solve_field; then
    echo "    solve-field found -- using the Astrometry.net binary."
    echo "    (not installing the PyPI 'astrometry' package: same import name)"
else
    # conda-forge builds `astrometry` for linux-64 and osx-64 ONLY. Asking for
    # it anywhere else fails the whole transaction, so the platform is checked
    # first rather than letting conda report a confusing solver error.
    case "${UNAME_S}/${UNAME_M}" in
        Linux/x86_64|Darwin/x86_64) CONDA_SOLVER_OK=1 ;;
        *)                          CONDA_SOLVER_OK=0 ;;
    esac

    if [[ "$MODE" != "venv" && $CONDA_SOLVER_OK -eq 1 ]] && command -v conda >/dev/null 2>&1; then
        echo "    Installing Astrometry.net from conda-forge."
        if [[ "$MODE" == "conda" ]]; then
            conda install -y -n "$ENV_NAME" -c conda-forge astrometry || true
        else
            conda install -y -c conda-forge astrometry || true
        fi
    fi

    if has_solve_field; then
        echo "    solve-field installed."
    else
        # The in-process solver covers every platform we support, including
        # Apple Silicon, where conda-forge has no astrometry build at all.
        case "${UNAME_S}/${UNAME_M}" in
            Linux/x86_64|Darwin/x86_64|Darwin/arm64)
                echo "    Using the in-process solver (no solve-field binary here)."
                "$PY" -m pip install -e ".[solver]"
                ;;
            *)
                warn "No plate solver is available for ${UNAME_S}/${UNAME_M}: conda-forge"
                warn "has no astrometry build and there is no in-process solver wheel."
                warn "Phase 2 will stack but cannot solve a WCS, so Phase 3 cannot"
                warn "measure a zero point. Run this inside x86-64 Linux (WSL works)."
                ;;
        esac
    fi
fi

# --- prove it works ----------------------------------------------------------
if [[ $RUN_DOCTOR -eq 1 ]]; then
    say "Checking the installation"
    "$PY" -m cassa_photometry.doctor || true
fi

say "Done"
case "$MODE" in
current) echo "    The pipeline is installed in '$CONDA_DEFAULT_ENV'. Try:  cassa --help" ;;
conda)   echo "    Activate it with:  conda activate $ENV_NAME"
         echo "    Then try:          cassa --help" ;;
venv)    echo "    Activate it with:  source .venv/bin/activate"
         echo "    Then try:          cassa --help" ;;
esac
