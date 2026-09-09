"""Supported platforms: Linux and macOS verified, Windows provisional.

The stance changed once already, so these tests pin it in the places a user
actually looks -- the doctor, the packaging metadata, the installer -- rather
than only in prose that can drift.

Windows was originally a hard failure because no plate solver was published for
it. ASTAP changed that: it ships command-line builds for win64, win32 and ARM64,
and every runtime dependency has a Windows wheel, so ``install.ps1`` can and
does install a working environment. What is still missing is evidence that it
works end to end, so the doctor reports ``WARN`` -- "this may work and nobody
has checked" -- rather than ``OK`` or ``FAIL``. Both of those would be a claim
the project cannot currently make.

WSL is real x86-64 Linux, so it is not a separate target: everything here
applies unchanged inside it, and it stays the recommended route until a native
install is confirmed.
"""

import re
from pathlib import Path
from unittest import mock

import pytest

from cassa_photometry import doctor

REPO = Path(__file__).resolve().parents[1]


# --- cassa-doctor states the position, before a run could hit it --------------

def test_native_windows_is_provisional_not_a_failure():
    """WARN, not FAIL: the solver blocker is gone, so refusing outright would
    be wrong -- but nothing has verified it, so passing silently would be too."""
    with mock.patch("platform.system", return_value="Windows"), \
         mock.patch("platform.machine", return_value="AMD64"):
        check = doctor.check_platform()

    assert check.status == doctor.WARN
    assert "provisional" in check.detail.lower()
    # It must name both routes: the one to try, and the one known to work.
    assert "install.ps1" in check.fix
    assert "WSL" in check.fix


def test_windows_does_not_make_the_doctor_exit_non_zero():
    """`cassa-doctor` exits on FAIL only. A Windows user with a working install
    must not be told their machine failed a check it did not fail."""
    with mock.patch("platform.system", return_value="Windows"), \
         mock.patch("platform.machine", return_value="AMD64"):
        assert doctor.check_platform().status != doctor.FAIL


@pytest.mark.parametrize(("system", "machine"), [
    ("Linux", "x86_64"),
    ("Linux", "aarch64"),
    ("Darwin", "arm64"),
])
def test_verified_platforms_pass(system, machine):
    with mock.patch("platform.system", return_value=system), \
         mock.patch("platform.machine", return_value=machine):
        assert doctor.check_platform().status == doctor.OK


def test_the_platform_check_runs_as_part_of_the_doctor():
    """A check nothing calls is not a check."""
    names = [c.name for c in doctor.run_checks(skip_network=True)]
    assert "platform" in names


def test_windows_solver_advice_does_not_point_at_apt():
    """The Linux fixes are actively misleading on Windows: `solve-field` is not
    published for it at all, so "not on PATH" would invite a hunt for something
    that does not exist."""
    with mock.patch("platform.system", return_value="Windows"):
        checks = {c.name: c for c in doctor.check_solvers()}

    assert "not published for Windows" in checks["solver: solve-field"].detail
    assert "apt" not in (checks["solver: astap"].fix or "")
    assert "install.ps1" in checks["solver: astap"].fix


# --- The packaging states it --------------------------------------------------

def test_packaging_declares_every_platform_the_installer_supports():
    text = (REPO / "pyproject.toml").read_text()
    classifiers = re.search(r"classifiers = \[(.*?)\]", text, re.S).group(1)

    assert "POSIX :: Linux" in classifiers
    assert "MacOS" in classifiers
    assert "Microsoft :: Windows" in classifiers, (
        "install.ps1 builds a working Windows environment, so the packaging "
        "should not claim otherwise"
    )


# --- install.ps1 actually installs --------------------------------------------

def test_install_ps1_installs_rather_than_redirecting():
    """It used to print WSL instructions and exit. It now has to do the work:
    an environment, the package, and the one solver Windows can run."""
    script = (REPO / "install.ps1").read_text()

    for required in ("pip install", "astap", "cassa_photometry.doctor"):
        assert required in script, f"install.ps1 never does {required!r}"


def test_install_ps1_still_offers_the_wsl_route():
    """WSL is the tested path, so it must stay one flag away."""
    script = (REPO / "install.ps1").read_text()
    assert "-Wsl" in script and "wsl --install" in script


def test_astap_download_asks_as_a_non_browser():
    """SourceForge serves a browser an HTML "your download will start shortly"
    page instead of the file, and Invoke-WebRequest's default User-Agent
    identifies as a browser -- which silently lands a ~114 KB HTML document
    named astap.zip. Measured 2026-09-09: default UA 114029 bytes of HTML,
    curl UA 322396 bytes beginning "PK"."""
    script = (REPO / "install.ps1").read_text()
    assert "-UserAgent" in script, (
        "without an explicit non-browser User-Agent the ASTAP download returns "
        "a web page, and every Windows install gets a broken solver"
    )


def test_the_astap_download_is_verified_before_it_is_trusted():
    """A guard against the host changing behaviour again: check the ZIP magic
    rather than the file extension, and fail loudly instead of leaving a broken
    install to surface at the first plate solve."""
    script = (REPO / "install.ps1").read_text()
    assert "0x50" in script and "0x4B" in script, (
        "install.ps1 should verify the download really is a ZIP (PK magic)"
    )


def test_the_installer_does_not_claim_windows_is_verified():
    """Listing Windows is honest; calling it tested is not, until it is."""
    script = (REPO / "install.ps1").read_text()
    assert "PROVISIONAL" in script or "provisional" in script
