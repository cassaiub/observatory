"""Raw-frame lookup for the phase 2 QA PDF.

The raw tree is a CLI argument and routinely lives outside the work directory,
so the QA PDF must not assume it sits beside the calibrated frames.
"""

import os

from astropy.io import fits

from cassa_photometry.phase2_integration.visuals import find_raw_frame


def _anchor(tmp_path, header=None):
    """A calibrated frame in <work>/phase1, as phase 1 would write it."""
    calibrated_dir = tmp_path / "work" / "phase1"
    calibrated_dir.mkdir(parents=True)
    path = calibrated_dir / "calibrated_frame_B.fits"
    fits.PrimaryHDU(data=None, header=header).writeto(path)
    return str(path)


def _raw(tmp_path, leaf="raw", name="frame_B.fits"):
    raw_dir = tmp_path / leaf
    raw_dir.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU().writeto(raw_dir / name)
    return str(raw_dir / name)


def test_header_provenance_finds_a_raw_tree_anywhere(tmp_path):
    """RAWFILE/RAWDIR win: the raw dir need not be near the work directory."""
    elsewhere = tmp_path / "archive" / "2023-08-23"
    elsewhere.mkdir(parents=True)
    fits.PrimaryHDU().writeto(elsewhere / "original.fits")

    header = fits.Header()
    header["RAWFILE"] = "original.fits"      # name need not match the calibrated one
    header["RAWDIR"] = str(elsewhere)
    anchor = _anchor(tmp_path, header)

    found, _ = find_raw_frame(anchor, os.path.basename(anchor), header)
    assert found == str(elsewhere / "original.fits")


def test_falls_back_to_conventional_layout_without_provenance(tmp_path):
    """Frames calibrated before RAWDIR existed still resolve via the layout."""
    expected = _raw(tmp_path, leaf="raw")     # sibling of <work>, lowercase
    anchor = _anchor(tmp_path)

    found, _ = find_raw_frame(anchor, os.path.basename(anchor), fits.Header())
    assert found == expected


def test_configured_raw_dir_overrides_stale_provenance(tmp_path):
    """A moved raw tree is recovered by phase2.raw_dir, which outranks RAWDIR."""
    expected = _raw(tmp_path, leaf="moved")
    header = fits.Header()
    header["RAWFILE"] = "frame_B.fits"
    header["RAWDIR"] = str(tmp_path / "gone")
    anchor = _anchor(tmp_path, header)

    found, _ = find_raw_frame(anchor, os.path.basename(anchor), header,
                              str(tmp_path / "moved"))
    assert found == expected


def test_missing_raw_reports_every_directory_it_tried(tmp_path):
    """A miss must name where it looked, or it cannot be acted on."""
    anchor = _anchor(tmp_path)
    found, searched = find_raw_frame(anchor, os.path.basename(anchor), fits.Header())
    assert found is None
    assert searched, "the warning would name no directories"
    assert len({os.path.normpath(d) for d in searched}) == len(searched)
