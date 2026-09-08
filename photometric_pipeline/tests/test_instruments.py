"""Instrument-profile registry, detector constants, and filter routing.

These cover the seam that makes the pipeline usable on more than one telescope:
the profile is selected by name, the header always outranks the profile's own
hardware knowledge, and every filter a setup can present routes to the right
reference band -- or to none at all, when no reference catalog describes it.
"""

import pytest
from astropy.io import fits
from conftest import FixtureProfile

from cassa_photometry.config import load_config
from cassa_photometry.instruments import (
    Cassa8InchProfile,
    InstrumentProfile,
    available_profiles,
    get_profile,
)


def _header(**cards):
    return fits.Header(cards)


# --- Registry -----------------------------------------------------------------

def test_registry_returns_the_named_profile():
    assert isinstance(get_profile("generic"), InstrumentProfile)
    assert isinstance(get_profile("cassa8"), Cassa8InchProfile)


def test_registry_ships_only_setups_this_observatory_operates():
    """Adding a telescope must not mean adding a line to this package."""
    assert available_profiles() == ["cassa8", "generic"]


def test_registry_defaults_to_generic():
    assert type(get_profile()) is InstrumentProfile
    assert type(get_profile(None)) is InstrumentProfile


def test_registry_is_case_insensitive():
    assert isinstance(get_profile("CASSA8"), Cassa8InchProfile)


def test_unknown_profile_names_the_alternatives():
    with pytest.raises(KeyError) as exc:
        get_profile("no_such_telescope")
    for name in available_profiles():
        assert name in str(exc.value)


def test_config_carries_the_profile_name():
    assert load_config().instrument == "generic"


def test_config_yaml_selects_a_profile(tmp_path):
    yaml_file = tmp_path / "cfg.yaml"
    yaml_file.write_text("instrument: cassa8\n")
    assert load_config(str(yaml_file)).instrument == "cassa8"


# --- The header outranks the profile ------------------------------------------

@pytest.mark.parametrize("name", ["generic", "cassa8"])
def test_header_wins_over_profile_knowledge(name):
    """A frame that states its own detector constants is believed, always."""
    header = _header(EGAIN=1.37, READNOIS=9.2, GAIN=100, SATURATE=61000)
    profile = get_profile(name)
    assert profile.get_gain(header) == pytest.approx(1.37)
    assert profile.get_read_noise(header) == pytest.approx(9.2)
    assert profile.get_saturation(header) == pytest.approx(61000)


def test_generic_admits_when_it_does_not_know():
    """No header, no hardware table: None, not a fabricated number."""
    profile = get_profile("generic")
    assert profile.get_gain(_header()) is None
    assert profile.get_read_noise(_header()) is None
    assert profile.get_saturation(_header()) is None


def test_unparseable_values_are_skipped_not_raised():
    """Cameras write text into numeric slots; fall through to the next key."""
    profile = get_profile("generic")
    assert profile.get_read_noise(_header(READNOIS="Mode0", RDNOISE=3.1)) == pytest.approx(3.1)


def test_a_profile_may_answer_where_the_generic_one_cannot():
    """The reason profiles exist: supplying what the header does not carry."""
    profile = FixtureProfile()
    assert profile.get_gain(_header()) == pytest.approx(1.0)
    assert profile.get_read_noise(_header()) == pytest.approx(10.0)
    # ...and it still yields to a header that states its own constants.
    assert profile.get_gain(_header(EGAIN=0.8)) == pytest.approx(0.8)


# --- CASSA 8-inch: gain-keyed CMOS detector curve -----------------------------

def test_cassa_reads_the_gain_curve_at_measured_points():
    profile = get_profile("cassa8")
    assert profile.get_gain(_header(GAIN=0)) == pytest.approx(0.58)
    assert profile.get_read_noise(_header(GAIN=100)) == pytest.approx(0.9)


def test_cassa_interpolates_between_measured_points():
    profile = get_profile("cassa8")
    read_noise = profile.get_read_noise(_header(GAIN=65))
    assert 0.9 < read_noise < 1.3


def test_cassa_resolves_the_hcg_transition():
    """Read noise drops sharply across the HCG switch; a scalar cannot say this."""
    profile = get_profile("cassa8")
    lcg = profile.get_read_noise(_header(GAIN=29))
    hcg = profile.get_read_noise(_header(GAIN=30))
    assert lcg > 2.0 and hcg < 1.5


def test_cassa_does_not_extrapolate_past_the_measured_range():
    profile = get_profile("cassa8")
    curve = profile.DETECTOR_CURVES["HCG"]
    assert profile.get_gain(_header(GAIN=10_000)) == pytest.approx(curve[max(curve)][0])


def test_cassa_never_reads_the_gain_setting_as_e_per_adu():
    """GAIN on this camera is the setting; misreading it is a ~100x error."""
    profile = get_profile("cassa8")
    assert profile.get_gain(_header(GAIN=100)) < 1.0


def test_cassa_knows_its_own_saturation():
    assert get_profile("cassa8").get_saturation(_header()) == pytest.approx(58000.0)


def test_cassa_plate_scale_matches_the_documented_optics():
    assert get_profile("cassa8").pixel_scale_arcsec() == pytest.approx(0.598, abs=0.001)


def test_detector_curve_is_replaceable_per_mode():
    """cassa-camchar output should not require editing source."""
    class Measured(Cassa8InchProfile):
        pass

    Measured.set_detector_curve({0: (0.5, 3.0), 100: (0.04, 0.8)}, mode="LCG")
    assert Measured.MEASURED is True
    assert Measured().get_gain(_header(GAIN=0, READOUTM="LCG")) == pytest.approx(0.5)
    # The other mode is left alone.
    assert Measured().get_gain(_header(GAIN=30, READOUTM="HCG")) == pytest.approx(0.14)
    # The shipped profile is untouched, and still flags itself as provisional.
    assert Cassa8InchProfile.MEASURED is False


# --- Readout mode -------------------------------------------------------------

def test_read_mode_is_taken_from_the_header_when_stated():
    profile = get_profile("cassa8")
    # Gain 100 would be inferred as HCG; an explicit LCG card overrules that.
    assert profile.read_mode(_header(GAIN=100, READOUTM="LCG")) == "LCG"
    assert profile.read_mode(_header(GAIN=0, READOUTM="HCG")) == "HCG"


def test_read_mode_is_inferred_from_gain_when_the_header_is_silent():
    profile = get_profile("cassa8")
    assert profile.read_mode(_header(GAIN=29)) == "LCG"
    assert profile.read_mode(_header(GAIN=30)) == "HCG"
    assert profile.read_mode(_header()) == "LCG"


@pytest.mark.parametrize("spelling,mode", [
    ("HCG", "HCG"), ("hcg", "HCG"), ("High Gain", "HCG"), ("HIGHGAIN", "HCG"),
    ("LCG", "LCG"), ("Low Conversion Gain", "LCG"), ("Standard", "LCG"),
])
def test_read_mode_spellings_are_normalised(spelling, mode):
    assert get_profile("cassa8").read_mode(_header(READOUTM=spelling)) == mode


def test_unknown_read_mode_falls_back_to_inference():
    """An unrecognised mode string must not crash or silently pick a curve."""
    profile = get_profile("cassa8")
    assert profile.read_mode(_header(GAIN=100, READOUTM="Mode7")) == "HCG"


def test_mode_changes_the_detector_constants_at_the_same_gain():
    """The whole point: one gain number, two different detectors."""
    profile = get_profile("cassa8")
    lcg = profile.get_read_noise(_header(GAIN=29, READOUTM="LCG"))
    hcg = profile.get_read_noise(_header(GAIN=29, READOUTM="HCG"))
    assert lcg > 2.0 and hcg < 1.5


def test_generic_profile_reports_the_raw_read_mode():
    assert get_profile("generic").get_read_mode(_header(READOUTM="HCG")) == "HCG"
    assert get_profile("generic").get_read_mode(_header()) is None


# --- Subframes ----------------------------------------------------------------

def test_subframe_origin_defaults_to_full_frame():
    assert get_profile("generic").get_subframe_origin(_header()) == (0, 0)


def test_subframe_origin_is_read_from_the_header():
    assert get_profile("cassa8").get_subframe_origin(
        _header(XORGSUBF=512, YORGSUBF=256)) == (512, 256)


def test_subframe_origin_accepts_the_alternate_spelling():
    assert get_profile("generic").get_subframe_origin(
        _header(SUBFRAMX=100, SUBFRAMY=200)) == (100, 200)


# --- Prior calibration --------------------------------------------------------

def test_raw_frame_is_not_flagged_as_calibrated():
    profile = get_profile("generic")
    assert profile.already_calibrated(_header(IMAGETYP="Light Frame")) is None
    assert profile.already_calibrated(_header(CALSTAT="")) is None
    assert profile.already_calibrated(_header(CALSTAT="None")) is None
    assert profile.already_calibrated(_header(BUNIT="ADU")) is None


@pytest.mark.parametrize("cards", [
    {"CALSTAT": "BDF"}, {"CALSTAT": "D"}, {"CALSTAT": "bdf"},
])
def test_calstat_flags_a_precalibrated_frame(cards):
    reason = get_profile("generic").already_calibrated(_header(**cards))
    assert reason is not None and "CALSTAT" in reason


def test_electron_units_flag_a_precalibrated_frame():
    """The pipeline's own output is in electrons; raw frames are in ADU."""
    reason = get_profile("generic").already_calibrated(_header(BUNIT="electron"))
    assert reason is not None and "BUNIT" in reason


# --- Filter routing -----------------------------------------------------------

@pytest.mark.parametrize("raw,band", [
    ("R", "R"), ("V", "V"), ("B", "B"), ("I", "I"), ("G", "G"),
    ("Red", "R"), ("Green", "G"), ("Blue", "B"),          # LRGB wheel
    ("Lum", "V"), ("L", "V"), ("Clear", "V"),
    ("V-Photometric", "V"), ("B-Photometric", "B"),       # long-form wheels
])
def test_science_bands_route_to_the_right_catalog(raw, band):
    assert get_profile("cassa8").science_band(_header(FILTER=raw)) == band


@pytest.mark.parametrize("raw", ["Ha", "SII", "OIII", "Dark", "DarkCap", "Blank"])
def test_uncalibratable_filters_report_no_band(raw):
    """No broadband catalog describes these, so no zero point may be invented."""
    assert get_profile("cassa8").science_band(_header(FILTER=raw)) is None


def test_missing_filter_falls_back_to_the_callers_default():
    """'The header did not say' is distinct from 'this filter has no band'."""
    assert get_profile("cassa8").science_band(_header(), default="V") == "V"


@pytest.mark.parametrize("raw,group", [
    ("R", "R_Photo"), ("Red", "Red"), ("Lum", "L"), ("Halpha", "Ha"), ("O-III", "OIII"),
])
def test_stacking_groups_are_canonical(raw, group):
    assert get_profile("cassa8").standardize_filter(_header(FILTER=raw)) == group


def test_bandwidth_suffix_separates_narrowband_stacks():
    profile = get_profile("cassa8")
    assert profile.standardize_filter(_header(FILTER="Ha", BANDWID=7)) == "Ha_7nm"


# --- Image classification -----------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("Bias Frame", "bias"), ("BIAS", "bias"),
    ("Dark Frame", "dark"), ("Flat Field", "flat"),
    ("Light Frame", "science"), ("", "science"),
])
def test_image_type_classification(raw, expected):
    assert get_profile("generic").get_image_type(_header(IMAGETYP=raw)) == expected


# --- Calibration conventions --------------------------------------------------

def test_flat_proxies_are_declared_by_the_profile_not_the_pipeline():
    assert get_profile("generic").flat_proxies() == {}
    assert get_profile("cassa8").flat_proxies() == {}
    # ...and a setup that does reuse one filter's flat declares it, rather than
    # the pipeline hardcoding any one observatory's habit.
    config = load_config()
    config.detector.flat_proxies = {"Red": "Luminance"}
    assert get_profile("generic", config=config).flat_proxies() == {"RED": "LUMINANCE"}
