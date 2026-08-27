import numpy as np

from cassa_photometry.phase3_photometry.engine import combine_zeropoints


def test_clean_zeropoints():
    zps = np.full(10, 25.0)
    errs = np.full(10, 0.1)
    zp, zp_err, n = combine_zeropoints(zps, errs)
    assert np.isclose(zp, 25.0)
    assert n == 10
    # Formal error = sqrt(1 / sum(1/err^2)) = sqrt(1 / (10 * 100)) ~= 0.0316
    assert np.isclose(zp_err, np.sqrt(1.0 / (10 * 100)), rtol=1e-3)


def test_outlier_is_rejected():
    zps = np.append(np.full(10, 25.0), 100.0)  # one gross outlier
    errs = np.full(11, 0.1)
    zp, zp_err, n = combine_zeropoints(zps, errs)
    assert n == 10                       # outlier clipped
    assert abs(zp - 25.0) < 0.01


def test_scatter_dominates_when_larger():
    # Widely scattered offsets but tiny formal errors -> scatter drives the error.
    zps = np.array([24.5, 25.0, 25.5])
    errs = np.array([0.001, 0.001, 0.001])
    _, zp_err, _ = combine_zeropoints(zps, errs, sigma=10.0)
    assert zp_err > 0.1
