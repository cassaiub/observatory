"""Supported platforms: Linux, macOS and Windows.

The stance has changed twice, so these tests pin it in the places a user
actually looks -- the doctor, the packaging metadata, the installer -- rather
than only in prose that can drift.

Windows was a hard failure for as long as no plate solver was published for it.
ASTAP ended that, and a native install has since been run on Windows: the
environment builds, the package imports, phase 1 reduces. Two Windows-only
defects surfaced in that run and are fixed -- see ``test_windows_file_handles``
and ``test_solver_index_use`` -- both invisible on POSIX, which is why they
lasted.

WSL remains a good route, and it is real x86-64 Linux so everything here applies
unchanged inside it. It is simply no longer the only one.
"""

import re
from pathlib import Path
from unittest import mock

import pytest

from cassa_photometry import doctor

REPO = Path(__file__).resolve().parents[1]


# --- cassa-doctor states the position, before a run could hit it --------------

@pytest.mark.parametrize(("system", "machine"), [
    ("Windows", "AMD64"),
    ("Windows", "ARM64"),
    ("Linux", "x86_64"),
    ("Linux", "aarch64"),
    ("Darwin", "arm64"),
])
def test_every_supported_platform_passes(system, machine):
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
    that does not exist.

    ASTAP's absence is pinned rather than inherited from the machine: on a
    developer box that has ASTAP the check is ``OK`` with no fix at all, and the
    assertion below would be testing nothing.
    """
    with mock.patch("platform.system", return_value="Windows"), \
         mock.patch("cassa_photometry.phase2_integration.solvers.astap.find_binary",
                    return_value=None):
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


def test_the_installer_registers_a_jupyter_kernel():
    """Installing packages sets up an environment; it does not make an existing
    Jupyter offer it. Without this the workshop notebooks fail with
    `ModuleNotFoundError: No module named 'cassa_photometry'`, which looks
    exactly like a failed install and is not one."""
    for script in ("install.sh", "install.ps1"):
        text = (REPO / script).read_text()
        assert "ipykernel install" in text, f"{script} never registers a kernel"
        assert "CASSA photometry" in text, f"{script} sets no kernel display name"
