import numpy as np
from astropy.io import fits

from cassa_photometry.fits_utils import (
    DQ_COSMIC_RAY,
    DQ_SATURATED,
    build_dq,
    read_mef,
    write_mef,
)


def test_build_dq_bits():
    sat = np.array([[True, False], [False, False]])
    cr = np.array([[True, False], [False, True]])
    dq = build_dq((2, 2), saturated=sat, cosmic_ray=cr)
    # Top-left has both flags OR'd together.
    assert dq[0, 0] == DQ_SATURATED | DQ_COSMIC_RAY
    assert dq[1, 1] == DQ_COSMIC_RAY
    assert dq[0, 1] == 0


def test_mef_roundtrip(tmp_path):
    sci = np.random.random((8, 8)).astype(np.float32)
    err = (np.random.random((8, 8)) + 0.1).astype(np.float32)
    dq = build_dq((8, 8), saturated=(sci > 0.9))
    header = fits.Header()
    header["OBJECT"] = "TESTFIELD"

    path = tmp_path / "out.fits"
    write_mef(str(path), sci=sci, err=err, dq=dq, header=header)

    r_sci, r_err, r_dq, r_hdr = read_mef(str(path))
    assert np.allclose(r_sci, sci, atol=1e-5)
    assert np.allclose(r_err, err, atol=1e-5)
    assert np.array_equal(r_dq, dq)
    assert r_hdr["OBJECT"] == "TESTFIELD"


def test_read_legacy_single_hdu(tmp_path):
    """A plain single-HDU FITS should read back with err/dq as None."""
    path = tmp_path / "legacy.fits"
    fits.PrimaryHDU(np.ones((4, 4), dtype=np.float32)).writeto(str(path))
    sci, err, dq, _ = read_mef(str(path))
    assert sci.shape == (4, 4)
    assert err is None and dq is None
