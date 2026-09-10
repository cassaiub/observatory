<#
.SYNOPSIS
    cassa-photometry installer for Windows.

.DESCRIPTION
    Picks an environment, installs every dependency including a working plate
    solver, installs the pipeline, and then runs `cassa-doctor` to prove it
    works. Running it twice is safe -- it updates in place. This is the
    PowerShell counterpart of install.sh, and does the same things in the same
    order.

    Native Windows is supported: ASTAP publishes command-line builds for
    win64, win32 and ARM64, and the runtime dependencies have Windows wheels.
    WSL also works if you already use it -- `-Wsl` prints those instructions
    instead.

.PARAMETER Conda
    Force the conda route.

.PARAMETER Venv
    Force a plain virtual environment in .\.venv

.PARAMETER NewEnv
    Build a dedicated environment even if one is already active.

.PARAMETER EnvName
    Conda environment name. Default: cassa-photometry

.PARAMETER NoDev
    Skip pytest and ruff.

.PARAMETER NoDoctor
    Skip the closing health check.

.PARAMETER Wsl
    Install nothing; print the WSL setup steps and exit.

.EXAMPLE
    .\install.ps1

.EXAMPLE
    .\install.ps1 -Venv -NoDev
#>

[CmdletBinding()]
param(
    [switch]$Conda,
    [switch]$Venv,
    [switch]$NewEnv,
    [string]$EnvName = "cassa-photometry",
    [switch]$NoDev,
    [switch]$NoDoctor,
    [switch]$Wsl
)

$ErrorActionPreference = "Stop"
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

function Head { param($m) Write-Host "`n==> $m" -ForegroundColor Cyan }
function Item { param($m) Write-Host "    $m" }
function Warn { param($m) Write-Host "[!] $m" -ForegroundColor Yellow }
function Die  { param($m) Write-Host "[x] $m" -ForegroundColor Red; exit 1 }

$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $RepoRoot

# --- the WSL route, on request -----------------------------------------------
if ($Wsl) {
    Head "cassa-photometry through WSL"
    Write-Host @"

  WSL is real x86-64 Linux, so every instruction in the documentation applies
  unchanged inside it. Native Windows works too -- run this script without
  -Wsl for that.
"@
    Head "1. Install WSL (once, from an admin PowerShell)"
    Item "wsl --install"
    Item ""
    Item "Reboot when it asks, then open 'Ubuntu' from the Start menu and let it"
    Item "finish creating your user account."

    Head "2. Install Miniforge inside WSL"
    Item 'wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"'
    Item "bash Miniforge3-Linux-x86_64.sh"

    Head "3. Clone and install, exactly as on Linux"
    Item "git clone https://github.com/cassaiub/observatory.git"
    Item "cd observatory/photometric_pipeline"
    Item "./install.sh"
    Item "conda activate cassa-photometry"
    Item "cassa-doctor"

    Head "Where your files live"
    Item "Your Windows drives are mounted under /mnt/c, /mnt/d and so on, so data"
    Item "you downloaded to Windows is reachable from WSL. Keep the repository and"
    Item "the work directory inside the Linux filesystem (your WSL home) though --"
    Item "reducing a night across /mnt is several times slower."
    exit 0
}

Head "cassa-photometry: installing natively on Windows"

# Wheels exist only for the Python versions each project has released for. On a
# newer Python, pip falls back to building from source -- and on Windows that
# needs the MSVC build tools, which most machines do not have. Saying so here
# costs nothing and saves a confusing compiler error later.
$pyMinor = [int](& { try { (python -c "import sys; print(sys.version_info.minor)") } catch { 0 } })
if ($pyMinor -ge 14) {
    Warn "This is Python 3.$pyMinor. Some dependencies (sep, photutils) publish no"
    Warn "wheel for it yet, so pip would build them from source -- which needs the"
    Warn "MSVC build tools. Two easier routes:"
    Item ""
    Item "  Use a Python that has wheels (3.10-3.13):"
    Item "    Remove-Item -Recurse -Force .venv; py -3.13 -m venv .venv; .\install.ps1 -Venv"
    Item ""
    Item "  Or use conda, which ships prebuilt binaries and needs no compiler:"
    Item "    .\install.ps1 -Conda"
    Item ""
}

# --- choose how to install ---------------------------------------------------
#
# Preference order, and why:
#   1. an already-active conda env  -- the user has told us where they work
#   2. conda, if available          -- prebuilt scientific stack, no compiler
#   3. a venv                       -- always works; every dependency has a wheel
$mode = ""
if ($Conda -and $Venv) { Die "-Conda and -Venv are mutually exclusive." }
if     ($Conda) { $mode = "conda" }
elseif ($Venv)  { $mode = "venv" }
else {
    $active = $env:CONDA_DEFAULT_ENV
    if ((-not $NewEnv) -and $active -and $active -ne "base") {
        $mode = "current"
    } elseif (Get-Command conda -ErrorAction SilentlyContinue) {
        $mode = "conda"
    } else {
        $mode = "venv"
    }
}

$py = ""

switch ($mode) {
    "current" {
        Head "Using the active conda environment: $($env:CONDA_DEFAULT_ENV)"
        Item "(pass -NewEnv to build a dedicated one instead)"
        $cmd = Get-Command python -ErrorAction SilentlyContinue
        if (-not $cmd) { Die "no python on PATH in the active environment." }
        $py = $cmd.Source
    }
    "conda" {
        if (-not (Get-Command conda -ErrorAction SilentlyContinue)) {
            Die "conda not found. Install Miniforge from https://conda-forge.org/download/, or re-run with -Venv."
        }
        $existing = (conda env list) | Select-String -Pattern "^\s*$([regex]::Escape($EnvName))\s"
        if ($existing) {
            Head "Updating conda environment '$EnvName'"
            conda env update -n $EnvName -f environment.yml --prune
        } else {
            Head "Creating conda environment '$EnvName'"
            conda env create -n $EnvName -f environment.yml
        }
        if ($LASTEXITCODE -ne 0) { Die "conda failed to build the environment." }
        $py = (conda run --no-capture-output -n $EnvName python -c "import sys; print(sys.executable)")
        $py = $py.Trim()
    }
    "venv" {
        $base = Get-Command python -ErrorAction SilentlyContinue
        if (-not $base) { $base = Get-Command py -ErrorAction SilentlyContinue }
        if (-not $base) { Die "no python found. Install Python 3.10+ from python.org, or Miniforge." }
        if (-not (Test-Path ".venv")) {
            Head "Creating virtual environment in .\.venv"
            & $base.Source -m venv .venv
        } else {
            Head "Reusing the virtual environment in .\.venv"
        }
        $py = Join-Path $RepoRoot ".venv\Scripts\python.exe"
        & $py -m pip install --quiet --upgrade pip setuptools wheel
    }
}

if (-not ($py -and (Test-Path $py))) { Die "could not locate the environment's python." }

# --- install the package -----------------------------------------------------
Head "Installing cassa-photometry (editable)"
$extras = if ($NoDev) { "" } else { "[dev]" }
& $py -m pip install -e ".$extras"
if ($LASTEXITCODE -ne 0) { Die "pip failed to install the package." }

# --- register a Jupyter kernel -----------------------------------------------
#
# Installing packages sets up an environment; it does not make Jupyter offer it.
# Someone who starts JupyterLab from an existing Anaconda, or opens a notebook in
# VS Code, sees only that Jupyter's kernels -- and the notebooks then fail with
# `ModuleNotFoundError: No module named 'cassa_photometry'`, which looks exactly
# like a failed install and is not one.
Head "Registering the Jupyter kernel"
& $py -c "import ipykernel" 2>$null
if ($LASTEXITCODE -eq 0) {
    & $py -m ipykernel install --user --name cassa-photometry `
          --display-name "Python (CASSA photometry)" 2>$null | Out-Null
    if ($LASTEXITCODE -eq 0) {
        Item "Kernel 'Python (CASSA photometry)' registered."
    } else {
        Warn "Could not register the kernel; the notebooks still run inside this env."
    }
} else {
    Item "ipykernel is not installed; skipping (only needed for the notebooks)."
    Item "Add it with:  $py -m pip install ipykernel"
}

# --- make sure there is a plate solver ---------------------------------------
#
# On Windows there is exactly one option, and that is why ASTAP exists in this
# pipeline at all: conda-forge builds astrometry.net for linux-64 and osx-64
# only, and the PyPI in-process solver publishes no Windows wheel. ASTAP does
# publish command-line builds for win64, win32 and ARM64, each of which is a
# single astap_cli.exe.
Head "Checking for a plate solver"
$envBin = Split-Path -Parent $py

function Find-Astap {
    $found = Get-Command astap_cli -ErrorAction SilentlyContinue
    if (-not $found) { $found = Get-Command astap -ErrorAction SilentlyContinue }
    if ($found) { return $found.Source }
    $local = Join-Path $envBin "astap_cli.exe"
    if (Test-Path $local) { return $local }
    return $null
}

function Get-AstapUrl {
    $base = "https://sourceforge.net/projects/astap-program/files/windows_installer"
    switch ($env:PROCESSOR_ARCHITECTURE) {
        "AMD64" { return "$base/astap_command-line_version_win64.zip/download" }
        "ARM64" { return "$base/astap_command-line_version_win11_aarch64.zip/download" }
        "x86"   { return "$base/astap_command-line_version_win32.zip/download" }
        default { return $null }
    }
}

function Install-Astap {
    # Returns $true only when astap_cli.exe is verifiably in place.
    #
    # Every command below is piped to Out-Null. A PowerShell function returns
    # everything written to its output stream, not just what `return` names, so
    # a single stray object would make the failure path truthy at the call site
    # -- reporting a solver that is not there. Expand-Archive in particular
    # emits an assembly reference the first time it runs.
    $url = Get-AstapUrl
    if (-not $url) {
        Warn "no ASTAP build is published for $($env:PROCESSOR_ARCHITECTURE)."
        return $false
    }

    $tmp = Join-Path ([IO.Path]::GetTempPath()) ("astap-" + [Guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $tmp -Force | Out-Null
    try {
        $zip = Join-Path $tmp "astap.zip"

        # The User-Agent is load-bearing. SourceForge serves a browser an HTML
        # "your download will start shortly" page instead of the file, and
        # Invoke-WebRequest's default User-Agent identifies as a browser -- so
        # the default lands a ~114 KB HTML document named astap.zip. Asking as
        # a non-browser client returns the archive itself.
        Invoke-WebRequest -Uri $url -OutFile $zip -UseBasicParsing `
                          -UserAgent "curl/8.0.0" | Out-Null

        # Verify it really is a ZIP rather than trusting the extension: if the
        # host ever changes behaviour again, fail here and say so, instead of
        # leaving a broken install to surface at the first plate solve.
        $magic = [IO.File]::ReadAllBytes($zip)[0..1]
        if ($magic[0] -ne 0x50 -or $magic[1] -ne 0x4B) {
            Warn "the download was not a ZIP archive (the host returned a web page)."
            Warn "Install ASTAP by hand from https://www.hnsky.org/astap.htm and set"
            Warn "phase2.astap_path to the astap_cli.exe you unpack."
            return $false
        }

        Expand-Archive -Path $zip -DestinationPath $tmp -Force | Out-Null
        $exe = @(Get-ChildItem -Path $tmp -Filter "astap*.exe" -Recurse) | Select-Object -First 1
        if (-not $exe) {
            Warn "the ASTAP archive held no executable."
            return $false
        }

        $dest = $envBin
        try {
            $probe = Join-Path $dest ".cassa-write-probe"
            [IO.File]::WriteAllText($probe, "x")
            Remove-Item $probe -Force
        } catch {
            $dest = Join-Path $env:LOCALAPPDATA "cassa-photometry\bin"
        }
        New-Item -ItemType Directory -Path $dest -Force | Out-Null

        $target = Join-Path $dest "astap_cli.exe"
        Copy-Item $exe.FullName $target -Force | Out-Null
        if (-not (Test-Path $target)) {
            Warn "could not place astap_cli.exe in $dest."
            return $false
        }

        $env:PATH = "$dest;$env:PATH"
        Item "ASTAP installed to $target"
        if ($dest -ne $envBin) {
            Warn "$dest is not permanently on your PATH. Either add it, or set"
            Warn "phase2.astap_path to $target in your config."
        }
        return $true
    } catch {
        Warn "could not download or unpack ASTAP: $($_.Exception.Message)"
        return $false
    } finally {
        Remove-Item $tmp -Recurse -Force -ErrorAction SilentlyContinue
    }
}

$astap = Find-Astap
if ($astap) {
    Item "ASTAP found at $astap -- using it."
} else {
    # Take the last value off the output stream and require it to be exactly
    # $true, so nothing incidental can be mistaken for success.
    $installed = @(Install-Astap)[-1] -eq $true
    if (-not $installed) {
        Warn "No plate solver could be installed."
        Warn "Phase 2 will stack but cannot solve a WCS, so Phase 3 cannot measure"
        Warn "a zero point. Install ASTAP by hand from https://www.hnsky.org/astap.htm"
        Warn "and set phase2.astap_path, or use WSL:  .\install.ps1 -Wsl"
    }
}

# --- prove it works ----------------------------------------------------------
if (-not $NoDoctor) {
    Head "Checking the installation"
    & $py -m cassa_photometry.doctor
}

Head "Done"
switch ($mode) {
    "current" { Item "The pipeline is installed in '$($env:CONDA_DEFAULT_ENV)'. Try:  cassa --help" }
    "conda"   { Item "Activate it with:  conda activate $EnvName"
                Item "Then try:          cassa --help" }
    "venv"    { Item "Activate it with:  .\.venv\Scripts\Activate.ps1"
                Item "Then try:          cassa --help" }
}
Write-Host ""
Item "If anything went wrong, 'cassa-doctor' prints the state of the whole"
Item "install in one screen -- that is the useful thing to send."
