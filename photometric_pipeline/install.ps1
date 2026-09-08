<#
.SYNOPSIS
    cassa-photometry on Windows: use WSL.

.DESCRIPTION
    There is no native Windows install, and this script does not attempt one.

    Phase 2 needs a plate solver, and neither channel ships one for Windows:
    conda-forge builds `astrometry` for linux-64 and osx-64 only, and the PyPI
    in-process solver publishes no Windows wheel. A native environment would
    therefore build, then fail at the first WCS solve -- and without a WCS,
    Phase 3 cannot cross-match a reference catalog, so there is no zero point
    and no calibrated magnitude.

    WSL gives you real x86-64 Linux, where both solvers work, and the Linux
    instructions then apply unchanged.

.EXAMPLE
    .\install.ps1
#>

$ErrorActionPreference = "Stop"

function Head { param($m) Write-Host "`n$m" -ForegroundColor Cyan }
function Item { param($m) Write-Host "  $m" }

Head "cassa-photometry: Windows setup goes through WSL"

Write-Host @"

  Phase 2 plate-solves, and no plate solver is published for native Windows --
  not on conda-forge, not on PyPI. Installing natively would give you an
  environment that builds and then cannot solve a WCS, which means no zero
  point and no calibrated magnitude.

  WSL is real x86-64 Linux, where the solver works and every other instruction
  applies unchanged.
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

Head "Where your files live"
Item "Your Windows drives are mounted under /mnt/c, /mnt/d and so on, so data"
Item "you downloaded to Windows is reachable from WSL. Keep the repository and"
Item "the work directory inside the Linux filesystem (your WSL home) though --"
Item "reducing a night across /mnt is several times slower."

Write-Host "`nFull instructions: workshop/handbook/participant_handbook.pdf`n" -ForegroundColor Cyan
