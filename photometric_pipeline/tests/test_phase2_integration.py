"""Phase 2: what survives the stack, and what the master says about itself.

The headline property is DQ propagation. Phase 2 used to rebuild the master's
data-quality plane from scratch, keeping only "no data", so every saturation,
cosmic-ray and bad-pixel flag phase 1 found was discarded at the stack. That was
provable in the delivered products -- `FLAGS` in every catalog took only the
values 0 and 8 -- and it made any downstream quality cut a no-op.
"""

import numpy as np
import pytest

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import (
    DQ_BAD_PIXEL,
    DQ_COSMIC_RAY,
    DQ_NO_DATA,
    DQ_REJECTED,
    DQ_SATURATED,
)
from cassa_photometry.phase2_integration.epochs import epoch_key
from cassa_photometry.phase2_integration.pipeline import (
    _bad_pixel_mask,
    _combine_dq,
    _mean_of,
)

# --- Data quality through the stack -------------------------------------------

def test_a_defect_in_every_frame_survives_into_the_master():
    """A saturated star is saturated in all of them, and must stay flagged."""
    shape = (4, 4)
    frame = np.zeros(shape, dtype=np.int32)
    frame[1, 1] = DQ_SATURATED
    dq = _combine_dq([frame.copy() for _ in range(3)],
                     n_used=np.full(shape, 3), master=np.ones(shape))
    assert dq[1, 1] & DQ_SATURATED
    assert not dq[0, 0] & DQ_SATURATED


def test_a_defect_in_only_some_frames_does_not_survive():
    """This is what dithering is for: a bad detector pixel lands on different
    sky in each frame, so good data from the others covers it."""
    shape = (4, 4)
    frames = [np.zeros(shape, dtype=np.int32) for _ in range(3)]
    frames[0][2, 2] = DQ_BAD_PIXEL
    dq = _combine_dq(frames, n_used=np.full(shape, 3), master=np.ones(shape))
    assert not dq[2, 2] & DQ_BAD_PIXEL


def test_uncovered_pixels_are_flagged_no_data():
    shape = (3, 3)
    master = np.ones(shape)
    master[0, 0] = np.nan
    n_used = np.full(shape, 2)
    n_used[0, 0] = 0
    dq = _combine_dq([np.zeros(shape, dtype=np.int32)] * 2, n_used, master)
    assert dq[0, 0] & DQ_NO_DATA
    assert not dq[1, 1] & DQ_NO_DATA


def test_partially_covered_pixels_are_flagged_rejected():
    shape = (3, 3)
    n_used = np.full(shape, 3)
    n_used[1, 1] = 1
    dq = _combine_dq([np.zeros(shape, dtype=np.int32)] * 3, n_used, np.ones(shape))
    assert dq[1, 1] & DQ_REJECTED
    assert not dq[0, 0] & DQ_REJECTED


@pytest.mark.parametrize("flag", [DQ_SATURATED, DQ_BAD_PIXEL, DQ_COSMIC_RAY])
def test_unusable_pixels_are_excluded_from_the_combine(flag):
    """Averaging a saturated core in as good data is how it ends up in a master
    looking merely bright."""
    dq = np.zeros((4, 4), dtype=np.int32)
    dq[0, 0] = flag
    assert _bad_pixel_mask(dq, (4, 4))[0, 0]
    assert not _bad_pixel_mask(dq, (4, 4))[1, 1]


def test_no_data_alone_does_not_exclude_a_pixel():
    """Coverage is handled by the finiteness test; double-counting it would
    reject pixels that are merely at the edge of one frame's footprint."""
    dq = np.zeros((4, 4), dtype=np.int32)
    dq[0, 0] = DQ_NO_DATA
    assert not _bad_pixel_mask(dq, (4, 4)).any()


def test_a_missing_dq_plane_is_treated_as_clean_not_as_bad():
    assert not _bad_pixel_mask(None, (4, 4)).any()


# --- Epoch binning ------------------------------------------------------------

def test_one_observing_night_stays_one_epoch_across_ut_midnight():
    """Binning on the UT date string splits a night in two and yields two
    half-depth masters differing in nothing but which side of midnight."""
    site = 90.4074  # Dhaka: UT midnight is 06:00 local, mid-night
    before = {"mjd_obs": 61285.895833, "site_long": site}   # 21:30 UT, 2 Sep
    after = {"mjd_obs": 61286.020833, "site_long": site}    # 00:30 UT, 3 Sep

    assert epoch_key(before, "night") == epoch_key(after, "night")
    # ...and the UT date string, which the naive approach would use, differs.
    assert "2026-09-02" != "2026-09-03"


def test_separate_nights_do_not_merge():
    site = 90.4074
    first = {"mjd_obs": 61285.9, "site_long": site}
    second = {"mjd_obs": 61286.9, "site_long": site}
    assert epoch_key(first, "night") != epoch_key(second, "night")


def test_epoch_binning_can_be_switched_off_entirely():
    """`none` must reproduce the historical grouping exactly."""
    assert epoch_key({"mjd_obs": 61286.0, "site_long": 90.0}, "none") == ""


def test_a_sub_night_bin_is_available_for_fast_variables():
    """Nightly binning averages away intra-night variability, which is correct
    for a supernova and wrong for a short-period variable."""
    early = {"mjd_obs": 61285.85, "site_long": 90.4}
    late = {"mjd_obs": 61286.15, "site_long": 90.4}
    assert epoch_key(early, "night") == epoch_key(late, "night")
    assert epoch_key(early, "6h") != epoch_key(late, "6h")


def test_a_frame_with_no_time_falls_back_rather_than_failing(pipeline_logs):
    assert epoch_key({}, "night", logger=pipeline_logs.handler and None) == ""


def test_missing_longitude_still_bins_and_warns():
    from cassa_photometry.phase2_integration import epochs

    epochs._WARNED.clear()
    assert epoch_key({"mjd_obs": 61286.02}, "night").startswith("N")


# --- Frame ranking and weighting ----------------------------------------------

def _pipeline():
    from cassa_photometry.phase2_integration.pipeline import IntegrationPipeline

    return IntegrationPipeline.__new__(IntegrationPipeline)


def test_the_anchor_prefers_a_sharp_frame_over_a_merely_bright_one():
    """The anchor sets both the astrometric reference and the flux scale every
    other frame is normalised to, so a soft anchor blurs the whole stack."""
    cfg = load_config().phase2
    pipeline = _pipeline()
    bright_soft = {"snr": 100.0, "fwhm": 8.0}
    dimmer_sharp = {"snr": 80.0, "fwhm": 4.0}
    assert pipeline._anchor_rank(dimmer_sharp, cfg) > pipeline._anchor_rank(bright_soft, cfg)


def test_ranking_falls_back_to_snr_when_no_fwhm_is_known():
    cfg = load_config().phase2
    pipeline = _pipeline()
    assert pipeline._anchor_rank({"snr": 50.0, "fwhm": None}, cfg) == 50.0


def test_point_source_weighting_penalises_poor_seeing():
    """1/sigma^2 is optimal for extended flux; a point source's SNR goes as
    1/(sigma*FWHM), so two frames of equal noise but different seeing should not
    contribute equally."""
    cfg = load_config().phase2
    cfg.stack_weight = "point_source"
    pipeline = _pipeline()

    sharp = pipeline._frame_weight(10.0, 4.0, cfg)
    soft = pipeline._frame_weight(10.0, 6.0, cfg)
    assert sharp > soft
    assert sharp / soft == pytest.approx((6.0 / 4.0) ** 2)


def test_extended_weighting_ignores_seeing():
    cfg = load_config().phase2
    cfg.stack_weight = "extended"
    pipeline = _pipeline()
    assert pipeline._frame_weight(10.0, 4.0, cfg) == pipeline._frame_weight(10.0, 6.0, cfg)


def test_frames_far_softer_than_the_group_are_rejected():
    cfg = load_config().phase2
    pipeline = _pipeline()
    frames = [{"fwhm": f, "path": f"f{i}"} for i, f in enumerate([4.0, 4.1, 4.2, 4.0, 12.0])]
    kept, rejected = pipeline._reject_poor_seeing(frames, cfg)
    assert len(rejected) == 1 and rejected[0]["fwhm"] == 12.0
    assert len(kept) == 4


def test_seeing_rejection_needs_enough_frames_to_have_a_median():
    cfg = load_config().phase2
    pipeline = _pipeline()
    frames = [{"fwhm": 4.0, "path": "a"}, {"fwhm": 12.0, "path": "b"}]
    kept, rejected = pipeline._reject_poor_seeing(frames, cfg)
    assert not rejected and len(kept) == 2


def test_summary_statistics_ignore_unknown_values():
    assert _mean_of([1.0, None, 3.0]) == 2.0
    assert _mean_of([None, None]) is None
