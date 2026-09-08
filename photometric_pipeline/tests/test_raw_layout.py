"""Reading the raw tree the acquisition software writes.

A night arrives sorted into ``<date>/BIAS|DARK|FLAT|LIGHT/<target>/``, not as
one flat directory of frames, and the pipeline is handed the top of that tree.
What sorts a frame is still its header -- the folders only say where files are.
"""

import glob
import os

import numpy as np
from astropy.io import fits

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import find_raw_frames, raw_tree_summary
from cassa_photometry.phase1_calibration.pipeline import run as run_phase1

LOGGER = get_logger("test_raw_layout")
SHAPE = (16, 16)


def _frame(directory, name, value, img_type, exptime, filt="R"):
    os.makedirs(directory, exist_ok=True)
    header = fits.Header({"IMAGETYP": img_type, "EXPTIME": exptime, "EGAIN": 1.0,
                          "READNOIS": 5.0, "FILTER": filt})
    path = os.path.join(directory, name)
    fits.PrimaryHDU(np.full(SHAPE, value, dtype=np.float32), header=header).writeto(path)
    return path


def _night(root, date="20260903", target="m22", n_science=2):
    """A night in the acquisition software's layout."""
    for i in range(2):
        _frame(os.path.join(root, date, "BIAS", "untargeted"),
               f"{date}T0545{i:02d}_untargeted_R_BIAS.fits", 500.0, "Bias Frame", 0.0)
    for i in range(n_science):
        _frame(os.path.join(root, date, "LIGHT", target),
               f"{date}T0437{i:02d}_{target}_R_LIGHT.fits", 1500.0, "Light Frame", 60.0)
    return root


def test_frames_are_found_below_the_directory_given(tmp_path):
    root = _night(str(tmp_path / "raw"))
    found = find_raw_frames(root)
    assert len(found) == 4
    assert found == sorted(found), "callers rely on a stable order"


def test_a_flat_directory_still_reads(tmp_path):
    """The simulator writes every frame side by side; that must keep working."""
    flat = tmp_path / "sim_raw"
    _frame(str(flat), "bias_001.fits", 500.0, "Bias", 0.0)
    _frame(str(flat), "sci_001.fits", 1500.0, "Light", 60.0)
    assert len(find_raw_frames(str(flat))) == 2


def test_a_single_file_is_returned_as_itself(tmp_path):
    path = _frame(str(tmp_path), "one.fits", 1500.0, "Light", 60.0)
    assert find_raw_frames(path) == [os.path.abspath(path)]


def test_the_summary_reports_each_directory(tmp_path):
    root = _night(str(tmp_path / "raw"))
    summary = raw_tree_summary(find_raw_frames(root), root)
    assert summary == {
        os.path.join("20260903", "BIAS", "untargeted"): 2,
        os.path.join("20260903", "LIGHT", "m22"): 2,
    }


def test_phase1_reduces_a_night_given_its_root(tmp_path):
    """The whole point: hand over the tree, get calibrated science frames."""
    root = _night(str(tmp_path / "raw"))
    out = str(tmp_path / "out")
    run_phase1(root, out, config=load_config(), logger=LOGGER)
    assert len(glob.glob(os.path.join(out, "calibrated_*.fits"))) == 2


def test_the_header_classifies_a_misfiled_frame_not_the_folder(tmp_path):
    """A bias filed under LIGHT/ is still a bias, because IMAGETYP says so."""
    root = str(tmp_path / "raw")
    _night(root)
    _frame(os.path.join(root, "20260903", "LIGHT", "m22"),
           "misfiled_bias.fits", 500.0, "Bias Frame", 0.0)

    out = str(tmp_path / "out")
    run_phase1(root, out, config=load_config(), logger=LOGGER)
    written = {os.path.basename(p) for p in glob.glob(os.path.join(out, "calibrated_*.fits"))}
    assert "calibrated_misfiled_bias.fits" not in written
    assert len(written) == 2


def test_repeated_filenames_across_nights_do_not_overwrite_each_other(tmp_path):
    """Two nights, same frame names: three frames in must be three frames out."""
    root = str(tmp_path / "raw")
    _night(root, date="20260903")
    _night(root, date="20260904")

    out = str(tmp_path / "out")
    run_phase1(root, out, config=load_config(), logger=LOGGER)
    written = glob.glob(os.path.join(out, "calibrated_*.fits"))
    assert len(written) == 4, "a name collision silently dropped a frame"


def test_provenance_records_the_subdirectory_a_frame_came_from(tmp_path):
    """RAWDIR must be the frame's own directory, not the root that was passed."""
    root = _night(str(tmp_path / "raw"))
    out = str(tmp_path / "out")
    run_phase1(root, out, config=load_config(), logger=LOGGER)

    calibrated = sorted(glob.glob(os.path.join(out, "calibrated_*.fits")))[0]
    header = fits.getheader(calibrated)
    assert os.path.isfile(os.path.join(header["RAWDIR"], header["RAWFILE"]))
    assert header["RAWDIR"].endswith(os.path.join("LIGHT", "m22"))
