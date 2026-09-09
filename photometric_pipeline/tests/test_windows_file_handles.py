"""Reading a frame must not keep its file open.

Astropy memory-maps image data by default, and a memory-mapped array keeps the
underlying file handle alive for as long as the array exists -- even after the
``HDUList`` has been closed. POSIX lets you unlink and replace an open file, so
on Linux and macOS this is invisible. Windows does not, and phase 2 is exactly
the shape that trips it: read a master stack, solve its WCS, write it back over
itself. Reported from a real Windows run as::

    PermissionError: [WinError 32] The process cannot access the file because
    it is being used by another process:
    '...\\work\\phase2\\Master_NGC7331_B_Photo_..._5fr.fits'

These tests encode the property rather than the platform, so they fail on any
machine if the policy regresses -- which matters, because the symptom only
appears on the one platform CI is least likely to run.
"""

import mmap

import numpy as np
import pytest
from astropy.io import fits

from cassa_photometry.fits_utils import open_fits, read_mef, write_mef

# Big enough that astropy actually memory-maps it. A tiny array can come back
# unmapped by chance, which would make these tests pass for the wrong reason.
SHAPE = (512, 512)


def _mmap_backed(array):
    """True when anything in the numpy base chain is a memory map."""
    node = array
    while node is not None:
        if isinstance(node, mmap.mmap):
            return True
        node = getattr(node, "base", None)
    return False


@pytest.fixture
def master(tmp_path):
    path = tmp_path / "Master_test.fits"
    write_mef(str(path),
              sci=np.ones(SHAPE, np.float32),
              err=np.full(SHAPE, 0.5, np.float32),
              dq=np.zeros(SHAPE, np.int32),
              header=fits.Header())
    return str(path)


def test_astropys_default_really_does_memory_map(master):
    """Guards the premise. If astropy stopped mapping by default these tests
    would pass while proving nothing, so assert the hazard still exists."""
    with fits.open(master) as hdul:
        data = hdul[0].data
    assert _mmap_backed(data), (
        "astropy no longer memory-maps by default; the tests below need "
        "rewriting rather than trusting"
    )


def test_open_fits_does_not_memory_map(master):
    with open_fits(master) as hdul:
        assert not _mmap_backed(hdul[0].data)
        assert not _mmap_backed(hdul["ERR"].data)
        assert not _mmap_backed(hdul["DQ"].data)


def test_read_mef_returns_arrays_that_own_their_memory(master):
    sci, err, dq, _ = read_mef(master)
    for name, plane in (("sci", sci), ("err", err), ("dq", dq)):
        assert not _mmap_backed(plane), f"{name} still references the file"


def test_a_frame_can_be_overwritten_immediately_after_being_read(master):
    """The phase 2 sequence, in miniature. On Windows this raises WinError 32
    if anything still holds the file."""
    sci, err, dq, header = read_mef(master)
    with open_fits(master) as hdul:
        solved = hdul[0].data

    write_mef(master, sci=solved, err=err, dq=dq, header=header,
              history="PHASE 2: WCS solved")

    again, _, _, _ = read_mef(master)
    assert again.shape == SHAPE


def test_every_read_in_the_package_goes_through_open_fits():
    """One policy, not a convention to remember. A bare ``fits.open`` would
    reintroduce the mapping silently, and only on Windows."""
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[1] / "src" / "cassa_photometry"
    offenders = []
    for path in root.rglob("*.py"):
        if path.name == "fits_utils.py":
            continue  # defines open_fits, and is allowed to call fits.open
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if "fits.open(" in line and not line.lstrip().startswith("#"):
                offenders.append(f"{path.relative_to(root)}:{number}")
    assert not offenders, (
        "use open_fits() instead of fits.open() at: " + ", ".join(offenders)
    )
