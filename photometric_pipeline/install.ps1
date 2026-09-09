<#
.SYNOPSIS
    cassa-photometry on Windows: use WSL.

.DESCRIPTION
    There is no native Windows install, and this script does not attempt one.

    This is about what is tested, not about what is theoretically possible.
    ASTAP -- the plate solver the pipeline prefers -- does publish a Windows
    build. What is missing is everything around it: install.sh is a shell
    script, the packaging classifiers and the 474-test suite target Linux and
    macOS, and no native Windows configuration has been verified end to end.

    WSL gives you real x86-64 Linux, which is the route that IS tested, and the
    Linux instructions then apply unchanged.

.EXAMPLE
    .\install.ps1
#>

$ErrorActionPreference = "Stop"

function Head { param($m) Write-Host "`n$m" -ForegroundColor Cyan }
function Item { param($m) Write-Host "  $m" }

Head "cassa-photometry: Windows setup goes through WSL"

Write-Host @"

  The pipeline is developed and tested on Linux and macOS, and its installer is
  a shell script. Rather than leave you with a half-working native environment,
  this script points you at WSL.

  WSL is real x86-64 Linux. Every instruction in the documentation applies
  unchanged inside it, and it is the configuration that is actually tested.
"@

Head "1. Install WSL (once, from an admin PowerShell)"
Item "wsl --install"
Item ""
Item "Reboot when it asks, then open 'Ubuntu' from the Start menu and let it"
Item "finish creating your user account."

Head "2. Install Miniforge inside WSL"
Item 'wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"'
Item "bash Miniforge3-Linux-x86_64.sh"
Item ""
Item "Accept the defaults, then close and reopen the Ubuntu window."

Head "3. Clone and install, exactly as on Linux"
Item "git clone https://github.com/cassaiub/observatory.git"
Item "cd observatory/photometric_pipeline"
Item "./install.sh"
Item "conda activate cassa-photometry"
Item "cassa-doctor"

Head "Which plate solver you will get"
Item "install.sh tries ASTAP first, then Astrometry.net's solve-field, then an"
Item "in-process Python solver, and uses the first that installs. Your results"
Item "do not depend on which one you end up with: the pipeline measures the"
Item "astrometric residual itself, so every backend produces the same header."
Item ""
Item "There is nothing to download in advance either -- the sky data a field"
Item "needs (about 6 MB for ASTAP) is fetched during the solve and cached."

Head "Where your files live"
Item "Your Windows drives are mounted under /mnt/c, /mnt/d and so on, so data"
Item "you downloaded to Windows is reachable from WSL. Keep the repository and"
Item "the work directory inside the Linux filesystem (your WSL home) though --"
Item "reducing a night across /mnt is several times slower."

Write-Host "`nFull instructions: workshop/handbook/participant_handbook.pdf`n" -ForegroundColor Cyan
