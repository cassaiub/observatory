"""Star/galaxy separation.

Classifying on ellipticity -- which is what this replaced -- confuses *shape*
with *concentration*: a face-on elliptical is round and passes as a star, a
trailed star is elongated and fails. That leaked galaxies into the zero point,
biasing every magnitude in the frame.
"""

import numpy as np
import pytest

from cassa_photometry.phase3_photometry.morphology import (
    AMBIGUOUS,
    EDGE,
    EXTENDED,
    SATURATED,
    STAR,
    classify,
)


def _population(n_stars=60, n_galaxies=25, seed=0, star_scatter=0.03):
    rng = np.random.default_rng(seed)
    difference = np.concatenate([
        rng.normal(0.0, star_scatter, n_stars),
        rng.normal(0.55, 0.15, n_galaxies),
    ])
    snr = np.concatenate([
        rng.uniform(30, 400, n_stars), rng.uniform(30, 300, n_galaxies)
    ])
    magnitude = np.concatenate([
        rng.uniform(12, 16, n_stars), rng.uniform(13, 16, n_galaxies)
    ])
    return difference, magnitude, snr, n_stars


def test_point_and_extended_sources_are_separated():
    difference, magnitude, snr, n_stars = _population()
    result = classify(difference, magnitude, snr)

    assert np.count_nonzero(result.classes[:n_stars] == STAR) >= int(0.95 * n_stars)
    assert np.count_nonzero(result.classes[n_stars:] == EXTENDED) >= int(0.9 * 25)


def test_the_threshold_is_calibrated_from_the_frames_own_locus():
    """No fixed cut: the locus moves with the seeing, and a cut that suits one
    night silently misclassifies another."""
    difference, magnitude, snr, _ = _population()
    tight = classify(difference, magnitude, snr)
    # The same population with a broader stellar locus, as poorer seeing gives.
    difference_soft, magnitude_soft, snr_soft, n_stars = _population(star_scatter=0.10)
    soft = classify(difference_soft, magnitude_soft, snr_soft)

    assert soft.locus_scatter > tight.locus_scatter * 2
    # ...and the stars are still recovered despite the wider locus.
    assert np.count_nonzero(soft.classes[:n_stars] == STAR) >= int(0.9 * n_stars)


def test_the_locus_survives_galaxies_outnumbering_stars():
    """Sigma-clipping finds the stars because they are the tight peak."""
    rng = np.random.default_rng(1)
    difference = np.concatenate([rng.normal(0.0, 0.03, 20), rng.normal(0.6, 0.2, 80)])
    snr = np.full(100, 100.0)
    magnitude = np.full(100, 14.0)

    result = classify(difference, magnitude, snr)
    assert result.locus_centre == pytest.approx(0.0, abs=0.05)
    assert np.count_nonzero(result.classes[:20] == STAR) >= 18


def test_a_sparse_field_returns_ambiguous_rather_than_a_silent_fixed_cut(pipeline_logs):
    """A high-latitude field may not define a locus at all. Guessing is worse
    than saying so."""
    result = classify(np.array([0.0, 0.1, 0.2]), np.array([14.0, 15.0, 16.0]),
                      np.array([100.0, 90.0, 80.0]), logger=pipeline_logs.handler and None)
    assert result.method == "insufficient-locus"
    assert set(result.classes) == {AMBIGUOUS}


def test_classification_reports_how_deep_it_can_be_trusted():
    """CLASSLIM, rather than pretending to classify to the detection limit."""
    difference, magnitude, snr, _ = _population()
    result = classify(difference, magnitude, snr)
    assert np.isfinite(result.limit_mag)
    assert 12.0 < result.limit_mag < 20.0


def test_sources_too_faint_to_classify_are_ambiguous():
    rng = np.random.default_rng(2)
    difference = np.concatenate([rng.normal(0.0, 0.03, 40), rng.normal(0.0, 0.4, 20)])
    snr = np.concatenate([rng.uniform(50, 400, 40), rng.uniform(3, 8, 20)])
    magnitude = np.concatenate([rng.uniform(12, 15, 40), rng.uniform(19, 21, 20)])

    result = classify(difference, magnitude, snr)
    assert np.count_nonzero(result.classes[40:] == AMBIGUOUS) >= 15


def test_saturated_and_edge_sources_override_any_measurement():
    """Their photometry is void, so a class derived from it would be fiction."""
    difference, magnitude, snr, _ = _population()
    saturated = np.zeros(difference.size, dtype=bool)
    saturated[0] = True
    edge = np.zeros(difference.size, dtype=bool)
    edge[1] = True

    result = classify(difference, magnitude, snr, saturated=saturated, edge=edge)
    assert result.classes[0] == SATURATED
    assert result.classes[1] == EDGE


def test_stellarity_is_continuous_and_bounded():
    difference, magnitude, snr, n_stars = _population()
    result = classify(difference, magnitude, snr)
    finite = np.isfinite(result.stellarity)

    assert np.all(result.stellarity[finite] >= 0.0)
    assert np.all(result.stellarity[finite] <= 1.0)
    assert np.nanmedian(result.stellarity[:n_stars]) > np.nanmedian(
        result.stellarity[n_stars:]
    )


def test_a_source_far_on_the_compact_side_is_not_called_a_star():
    """More compact than the PSF is not a better star -- it is an artefact."""
    difference, magnitude, snr, _ = _population()
    difference = np.append(difference, -1.5)      # impossibly compact
    magnitude = np.append(magnitude, 14.0)
    snr = np.append(snr, 200.0)

    result = classify(difference, magnitude, snr)
    assert result.classes[-1] == AMBIGUOUS


def test_sources_with_no_measurement_are_ambiguous():
    difference, magnitude, snr, _ = _population()
    difference = np.append(difference, np.nan)
    magnitude = np.append(magnitude, 15.0)
    snr = np.append(snr, 50.0)

    result = classify(difference, magnitude, snr)
    assert result.classes[-1] == AMBIGUOUS
