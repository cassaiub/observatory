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
# Phase 2 needs a plate solver, and three backends can provide one. They are
# tried in this order, and the first that installs wins:
#
#   1. ASTAP        -- a sub-megabyte binary with no Python dependency, packaged
#                      for Debian/Ubuntu and published for Linux x86-64 and
#                      aarch64 and for macOS on Intel and Apple Silicon. It is
#                      the only backend that exists on ARM Linux at all, and it
#                      fetches just the sky tiles a field needs (~6 MB) instead
#                      of a multi-gigabyte index set.
#   2. solve-field  -- the Astrometry.net binary, from conda-forge. The
#                      reference implementation, and the only backend that
#                      supplies its own matched-star table.
#   3. astrometry   -- the PyPI in-process solver, so a pip-only machine with no
#                      system packages still works.
#
# solve-field and the PyPI package MUST NOT coexist: conda-forge's astrometry
# package installs Python bindings that import under the same name, and
# whichever lands second wins. Step 3 therefore only runs when step 2 did not.
#
# Whichever lands, the pipeline produces the same products: the astrometric
# residual is measured against a reference catalog when the backend cannot
# supply matched stars, so ASTRMS exists on every path.
say "Checking for a plate solver"
ENV_BIN="$(dirname "$PY")"
UNAME_S="$(uname -s 2>/dev/null || echo unknown)"
UNAME_M="$(uname -m 2>/dev/null || echo unknown)"

has_astap() {
    command -v astap_cli >/dev/null 2>&1 || command -v astap >/dev/null 2>&1
}

has_solve_field() {
    command -v "$ENV_BIN/solve-field" >/dev/null 2>&1 || command -v solve-field >/dev/null 2>&1
}

# The upstream command-line zip, per platform. Empty means "no build we can
# fetch unattended"; that platform falls through to the next backend.
astap_zip_url() {
    local base="https://sourceforge.net/projects/astap-program/files"
    case "${UNAME_S}/${UNAME_M}" in
        Linux/x86_64)   echo "$base/linux_installer/astap_command-line_version_Linux_amd64.zip/download" ;;
        Linux/aarch64)  echo "$base/linux_installer/astap_command-line_version_Linux_aarch64.zip/download" ;;
        Darwin/x86_64)  echo "$base/macOS%20installer/astap_command-line_version_macOS_x86_64.zip/download" ;;
        Darwin/arm64)   echo "$base/macOS%20installer/astap_command-line_version_macOS_M1.zip/download" ;;
        *)              echo "" ;;
    esac
}

install_astap() {
    # Prefer the distribution package: it is signed, updated with the system,
    # and needs no unpacking. Only fall back to the upstream zip when apt is
    # absent or the package is not there.
    if command -v apt-get >/dev/null 2>&1 && [[ $(id -u) -eq 0 || -n "${SUDO_USER:-}" ]]; then
        if apt-get install -y astap-cli >/dev/null 2>&1; then
            has_astap && { echo "    ASTAP installed from apt (astap-cli)."; return 0; }
        fi
    fi

    local url
    url="$(astap_zip_url)"
    [[ -n "$url" ]] || return 1
    command -v curl >/dev/null 2>&1 || return 1
    command -v unzip >/dev/null 2>&1 || return 1

    local dest="$ENV_BIN"
    [[ -w "$dest" ]] || dest="$HOME/.local/bin"
    mkdir -p "$dest" || return 1

    local tmp
    tmp="$(mktemp -d)" || return 1
    if curl -fsSL -o "$tmp/astap.zip" "$url" && unzip -o -q "$tmp/astap.zip" -d "$tmp"; then
        # The zip holds the bare executable; name varies slightly by platform.
        local binary
        binary="$(find "$tmp" -maxdepth 2 -type f -name 'astap*' ! -name '*.zip' | head -1)"
        if [[ -n "$binary" ]]; then
            install -m 0755 "$binary" "$dest/astap_cli" 2>/dev/null || {
                cp "$binary" "$dest/astap_cli" && chmod +x "$dest/astap_cli"; }
            rm -rf "$tmp"
            export PATH="$dest:$PATH"
            has_astap && { echo "    ASTAP installed to $dest/astap_cli."; return 0; }
        fi
    fi
    rm -rf "$tmp"
    return 1
}

if has_astap; then
    echo "    ASTAP found -- using it."
elif install_astap; then
    :
else
    echo "    ASTAP unavailable here; trying Astrometry.net."

    if ! has_solve_field; then
        # conda-forge builds `astrometry` for linux-64 and osx-64 ONLY. Asking
        # for it anywhere else fails the whole transaction, so the platform is
        # checked first rather than letting conda report a confusing error.
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
    fi

    if has_solve_field; then
        echo "    solve-field available."
    else
        # The in-process solver covers Linux x86-64 and both macOS
        # architectures. It publishes no wheel for ARM Linux, which is exactly
        # where ASTAP is the only option.
        case "${UNAME_S}/${UNAME_M}" in
            Linux/x86_64|Darwin/x86_64|Darwin/arm64)
                echo "    Using the in-process solver (no ASTAP or solve-field here)."
                "$PY" -m pip install -e ".[solver]"
                ;;
            *)
                warn "No plate solver could be installed for ${UNAME_S}/${UNAME_M}."
                warn "Phase 2 will stack but cannot solve a WCS, so Phase 3 cannot"
                warn "measure a zero point. Install ASTAP by hand from"
                warn "https://www.hnsky.org/astap.htm and set phase2.astap_path."
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
