import numpy as np
import pytest
from astropy.io import fits
from conftest import FixtureProfile

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
    std_ccd = load_standardized_ccds(path, FixtureProfile())[0]
    expected = np.sqrt(2.0 * 1000.0 + 5.0 ** 2) / 2.0  # = 22.5
    assert std_ccd.ccd.uncertainty is not None
    assert np.allclose(std_ccd.ccd.uncertainty.array, expected, rtol=1e-6)
    assert not std_ccd.sat_mask.any()  # nothing near saturation


def test_weighted_stack_propagates_variance():
    """Weighted mean with propagated variance: sum(w^2 var)/(sum w)^2."""
    cube = np.array([[[10.0]], [[20.0]]])   # two 1x1 frames
    varcube = np.array([[[4.0]], [[4.0]]])
    weights = np.array([0.5, 0.5])
    master, master_var, n_used = MathEngine.weighted_stack(
        cube, varcube, weights, use_sigma_clip=True, sigma=5.0,
    )
    assert np.isclose(master[0, 0], 15.0)
    assert np.isclose(master_var[0, 0], 2.0)
    assert n_used[0, 0] == 2


def test_a_small_stack_keeps_every_frame():
    """Min/max rejection discarded two frames of every three: a 3-frame stack
    kept ONE frame per pixel, and where values tied it could reject them all
    and leave a NaN."""
    rng = np.random.default_rng(0)
    cube = rng.normal(100.0, 5.0, (3, 8, 8))
    varcube = np.full_like(cube, 25.0)

    master, master_var, n_used = MathEngine.weighted_stack(
        cube, varcube, np.ones(3) / 3, use_sigma_clip=False,
    )
    assert np.all(n_used == 3)
    assert np.isfinite(master).all() and np.isfinite(master_var).all()


def test_identical_frames_do_not_all_get_rejected():
    """The pathological case of the old min/max rule: every frame equals both
    the min and the max, so all were rejected and the pixel became NaN."""
    cube = np.full((3, 4, 4), 100.0)
    varcube = np.full((3, 4, 4), 25.0)

    master, _, n_used = MathEngine.weighted_stack(
        cube, varcube, np.ones(3) / 3, use_sigma_clip=False,
    )
    assert np.isfinite(master).all()
    assert np.all(n_used == 3)


def test_flagged_pixels_are_excluded_from_the_combine():
    cube = np.stack([np.full((4, 4), 100.0), np.full((4, 4), 100.0),
                     np.full((4, 4), 5000.0)])       # one saturated frame
    varcube = np.full_like(cube, 25.0)
    bad = np.zeros(cube.shape, dtype=bool)
    bad[2] = True

    master, _, n_used = MathEngine.weighted_stack(
        cube, varcube, np.ones(3) / 3, use_sigma_clip=False, bad_pixels=bad,
    )
    assert np.allclose(master, 100.0)
    assert np.all(n_used == 2)


def test_a_pixel_with_no_coverage_is_reported_not_invented():
    cube = np.full((2, 3, 3), 10.0)
    cube[:, 1, 1] = np.nan
    varcube = np.full_like(cube, 1.0)

    master, _, n_used = MathEngine.weighted_stack(cube, varcube, np.ones(2) / 2)
    assert n_used[1, 1] == 0
    assert np.isnan(master[1, 1])


def test_the_flux_scale_tracks_transparency_not_seeing():
    """A peak-pixel ratio reads a seeing change as a transparency change, and
    that spurious factor then multiplies the whole frame before stacking."""

    positions = np.array([[64.0, 64.0], [128.0, 100.0], [90.0, 150.0], [160.0, 60.0]])

    def frame(fwhm, throughput):
        image = np.zeros((200, 200))
        for x, y in positions:
            from cassa_photometry.simulate.scene import _add_moffat
            _add_moffat(image, x, y, 10000.0 * throughput, fwhm)
        return image

    reference = frame(4.0, 1.0)

    # The aperture is sized from the worse of the two frames' seeing.
    radius = 6.0 * 3.0

    # Same transparency, much worse seeing: the scale must stay ~1. The peak
    # pixel gives 2.25 here -- exactly (6/4)**2 -- so the frame would have been
    # multiplied by 2.25 before stacking, purely because the seeing worsened.
    same_flux = MathEngine.calc_scale(frame(6.0, 1.0), reference, positions, positions,
                                      aperture_radius=radius)
    assert same_flux == pytest.approx(1.0, abs=0.03)

    # Genuinely fainter frame at the same seeing: the scale must find it.
    dimmer = MathEngine.calc_scale(frame(4.0, 0.5), reference, positions, positions,
                                   aperture_radius=radius)
    assert dimmer == pytest.approx(2.0, rel=0.03)
