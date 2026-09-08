"""Supported platforms: Linux and macOS. Native Windows is not one of them.

This is a deliberate scope decision, not an oversight. Phase 2 plate-solves, and
neither channel publishes a solver for Windows -- conda-forge builds
``astrometry`` for linux-64 and osx-64 only, and the PyPI in-process solver ships
no Windows wheel. A native Windows install would therefore succeed and then fail
at the first WCS solve, with no zero point and no calibrated magnitude to show
for it.

WSL is real x86-64 Linux, so it is not a separate target: everything here
applies unchanged inside it. These tests exist so the stance stays stated in the
places a user actually looks -- the doctor, the packaging metadata, the
installer -- rather than only in prose that can drift.
"""

import re
from pathlib import Path
from unittest import mock

import pytest

from cassa_photometry import doctor

REPO = Path(__file__).resolve().parents[1]


# --- cassa-doctor says so, before the solver failure would ---------------------

def test_native_windows_is_reported_as_a_failure():
    with mock.patch("platform.system", return_value="Windows"), \
         mock.patch("platform.machine", return_value="AMD64"):
        check = doctor.check_platform()

    assert check.status == doctor.FAIL
    assert "not supported" in check.detail.lower()
    # The fix must name the actual route, not merely refuse.
    assert "WSL" in check.fix


@pytest.mark.parametrize(("system", "machine"), [
    ("Linux", "x86_64"),
    ("Darwin", "arm64"),
])
def test_supported_platforms_pass(system, machine):
    with mock.patch("platform.system", return_value=system), \
         mock.patch("platform.machine", return_value=machine):
        assert doctor.check_platform().status == doctor.OK


def test_the_platform_check_runs_as_part_of_the_doctor():
    """A check nothing calls is not a check."""
    names = [c.name for c in doctor.run_checks(skip_network=True)]
    assert "platform" in names


# --- The packaging states it ---------------------------------------------------

def test_packaging_declares_only_linux_and_macos():
    text = (REPO / "pyproject.toml").read_text()
    classifiers = re.search(r"classifiers = \[(.*?)\]", text, re.S).group(1)

    assert "POSIX :: Linux" in classifiers
    assert "MacOS" in classifiers
    assert "Microsoft" not in classifiers, (
        "a Windows classifier would advertise a configuration that cannot "
        "plate-solve"
    )


# --- The Windows entry point points at WSL and installs nothing ----------------

def test_install_ps1_installs_nothing_and_redirects_to_wsl():
    script = (REPO / "install.ps1").read_text()

    assert "WSL" in script
    # It must not try to build a native environment: no conda/pip/env creation.
    for forbidden in ("conda create", "conda env create", "pip install",
                      "python -m venv"):
        assert forbidden not in script, (
            f"install.ps1 attempts {forbidden!r}; native Windows is unsupported "
            f"and the script exists only to redirect to WSL"
        )


def test_the_installer_does_not_claim_windows_support():
    text = (REPO / "install.sh").read_text()
    assert "needed on Windows" not in text
