import numpy as np
from astropy.io import fits

from cassa_photometry.instruments import ITelescopeNetworkProfile
from cassa_photometry.phase1_calibration.data_models import load_standardized_ccds
from cassa_photometry.phase2_integration.math_utils import MathEngine


def _make_frame(tmp_path, value=1000.0, gain=2.0, readnoise=5.0):
    header = fits.Header()
    header["IMAGETYP"] = "Light"
    header["EXPTIME"] = 10.0
    header["EGAIN"] = gain
    header["READNOIS"] = readnoise
    header["FILTER"] = "R"
    header["INSTRUME"] = "T32"
    header["TELESCOP"] = "T32"
    data = np.full((16, 16), value, dtype=np.float32)
    path = tmp_path / "frame.fits"
    fits.PrimaryHDU(data=data, header=header).writeto(str(path))
    return str(path)


def test_create_deviation_seeds_error_budget(tmp_path):
    """sigma_ADU = sqrt(gain*ADU + readnoise^2) / gain."""
    path = _make_frame(tmp_path, value=1000.0, gain=2.0, readnoise=5.0)
    std_ccd = load_standardized_ccds(path, ITelescopeNetworkProfile())[0]
    expected = np.sqrt(2.0 * 1000.0 + 5.0 ** 2) / 2.0  # = 22.5
    assert std_ccd.ccd.uncertainty is not None
    assert np.allclose(std_ccd.ccd.uncertainty.array, expected, rtol=1e-6)
    assert not std_ccd.sat_mask.any()  # nothing near saturation


def test_weighted_stack_propagates_variance():
    """Weighted mean with propagated variance: sum(w^2 var)/(sum w)^2."""
    cube = np.array([[[10.0]], [[20.0]]])   # two 1x1 frames
    varcube = np.array([[[4.0]], [[4.0]]])
    weights = np.array([0.5, 0.5])
    master, master_var = MathEngine.weighted_stack(
        cube, varcube, weights, use_sigma_clip=True, sigma=5.0,
    )
    assert np.isclose(master[0, 0], 15.0)
    assert np.isclose(master_var[0, 0], 2.0)
