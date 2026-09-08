"""Describing an observing setup the pipeline ships no profile for.

The pipeline deliberately registers only the setups CASSA operates, so these
cover the three routes a new setup takes instead -- a config block, a local
profile file, an installed plugin -- and the plate-scale resolution that phase 2
and phase 3 depend on.
"""

import numpy as np
import pytest
from astropy.io import fits
from conftest import header

from cassa_photometry.config import ConfigError, load_config
from cassa_photometry.instruments import Cassa8InchProfile, get_profile

EXAMPLES = "examples/profiles"


def _config(**detector):
    config = load_config()
    for key, value in detector.items():
        setattr(config.detector, key, value)
    return config


# --- Route 1: a detector: block, no code --------------------------------------

def test_detector_block_supplies_what_the_header_lacks():
    config = _config(gain=1.4, read_noise=7.0, saturation_adu=60000, pixel_scale_arcsec=0.40)
    profile = get_profile("generic", config=config)

    assert profile.get_gain(header()) == pytest.approx(1.4)
    assert profile.get_read_noise(header()) == pytest.approx(7.0)
    assert profile.get_saturation(header()) == pytest.approx(60000)
    assert profile.get_pixel_scale(header()) == pytest.approx(0.40)


def test_the_header_still_outranks_a_detector_block():
    """The frame is the first authority on its own data."""
    config = _config(gain=1.4, read_noise=7.0)
    profile = get_profile("generic", config=config)
    frame = header(EGAIN=0.8, READNOIS=2.0)

    assert profile.get_gain(frame) == pytest.approx(0.8)
    assert profile.get_read_noise(frame) == pytest.approx(2.0)


def test_override_header_inverts_that_and_says_so(pipeline_logs):
    """Discarding a measured card for a typed one must never be silent."""
    from cassa_photometry.instruments import overrides

    overrides._WARNED.clear()
    config = _config(gain=1.4, override_header=True)
    profile = get_profile("generic", config=config)

    assert profile.get_gain(header(EGAIN=0.8)) == pytest.approx(1.4)
    assert "override_header" in pipeline_logs.text
    assert "0.8" in pipeline_logs.text


def test_filter_tables_are_merged_not_replaced():
    """Adding one odd filter name must not lose the standard ones."""
    config = _config(filter_map={"Sloan-R": "R_Photo"}, science_bands={"Sloan-R": "R"})
    profile = get_profile("generic", config=config)

    assert profile.standardize_filter(header(FILTER="Sloan-R")) == "R_Photo"
    assert profile.science_band(header(FILTER="Sloan-R")) == "R"
    assert profile.standardize_filter(header(FILTER="V")) == "V_Photo"


def test_uncalibrated_filters_can_be_extended():
    config = _config(uncalibrated_filters=["Custom-NB"])
    profile = get_profile("generic", config=config)

    assert profile.science_band(header(FILTER="Custom-NB")) is None
    assert profile.science_band(header(FILTER="HA")) is None  # still excluded


def test_a_detector_block_layers_over_any_profile():
    config = _config(read_noise=99.0)
    profile = get_profile("cassa8", config=config)

    assert isinstance(profile, Cassa8InchProfile)
    assert profile.get_read_noise(header()) == pytest.approx(99.0)
    # ...and the profile's own knowledge is untouched where not overridden.
    assert profile.get_saturation(header()) == pytest.approx(58000.0)


def test_no_detector_block_leaves_the_profile_alone():
    assert type(get_profile("cassa8", config=load_config())) is Cassa8InchProfile


def test_a_mistyped_detector_value_is_rejected(tmp_path):
    path = tmp_path / "cfg.yaml"
    path.write_text("detector:\n  gain: not-a-number\n")
    with pytest.raises(ConfigError) as exc:
        load_config(str(path))
    assert "detector.gain" in str(exc.value)


# --- Route 2: a profile class in a local file ---------------------------------

def test_instrument_module_loads_a_profile_from_a_path(tmp_path):
    path = tmp_path / "my_profile.py"
    path.write_text(
        "from cassa_photometry.instruments.base import InstrumentProfile\n"
        "class MyProfile(InstrumentProfile):\n"
        "    @property\n"
        "    def name(self):\n"
        "        return 'My Rig'\n"
        "    def get_gain(self, header):\n"
        "        return super().get_gain(header) or 2.5\n"
    )
    config = load_config()
    config.instrument_module = f"{path}:MyProfile"
    profile = get_profile(config.instrument, config=config)

    assert profile.name == "My Rig"
    assert profile.get_gain(header()) == pytest.approx(2.5)
    assert profile.get_gain(header(EGAIN=1.1)) == pytest.approx(1.1)


@pytest.mark.parametrize(
    "spec, expected_name",
    [
        (f"{EXAMPLES}/template.py:MyObservatoryProfile", "My Observatory"),
        (f"{EXAMPLES}/itelescope.py:ITelescopeNetworkProfile", "iTelescope"),
    ],
)
def test_the_shipped_examples_actually_load(spec, expected_name):
    """A template nobody can run is documentation, not an example."""
    config = load_config()
    config.instrument_module = spec
    assert expected_name in get_profile(config.instrument, config=config).name


def test_a_bad_instrument_module_says_what_is_wrong(tmp_path):
    config = load_config()

    config.instrument_module = "nowhere.py:Missing"
    with pytest.raises(KeyError, match="not found"):
        get_profile(None, config=config)

    path = tmp_path / "empty.py"
    path.write_text("x = 1\n")
    config.instrument_module = f"{path}:Absent"
    with pytest.raises(KeyError, match="no class named"):
        get_profile(None, config=config)

    config.instrument_module = str(path)
    with pytest.raises(KeyError, match="ClassName"):
        get_profile(None, config=config)


# --- Route 3: an installed plugin ---------------------------------------------

def test_entry_point_plugins_join_the_registry(monkeypatch):
    from cassa_photometry.instruments import base, registry

    class PluginProfile(base.InstrumentProfile):
        @property
        def name(self):
            return "Plugin Scope"

    monkeypatch.setattr(registry, "_plugin_profiles", lambda: {"plugin": PluginProfile})
    assert "plugin" in registry.available_profiles()
    assert registry.get_profile("plugin").name == "Plugin Scope"


def test_a_broken_plugin_does_not_break_the_pipeline(monkeypatch):
    """A third party's bug must not stop a reduction with a built-in profile."""
    class Boom:
        name = "broken"

        def load(self):
            raise ImportError("no")

    from cassa_photometry.instruments import registry

    monkeypatch.setattr(registry, "entry_points", None, raising=False)
    monkeypatch.setattr(
        "importlib.metadata.entry_points", lambda **kw: [Boom()]
    )
    assert registry.get_profile("generic") is not None


# --- Plate scale --------------------------------------------------------------

def test_pixel_scale_prefers_the_header():
    """SECPIX sits unread in many headers; reading it is what gives phase 2 a hint."""
    assert get_profile("generic").get_pixel_scale(header(SECPIX=0.4)) == pytest.approx(0.4)
    # The header wins even over a profile that knows its own optics.
    assert get_profile("cassa8").get_pixel_scale(header(SECPIX=0.4)) == pytest.approx(0.4)


def test_pixel_scale_from_pixel_size_and_focal_length():
    scale = get_profile("generic").get_pixel_scale(header(XPIXSZ=2.9, FOCALLEN=1000.0))
    assert scale == pytest.approx(2.9 / 1000.0 * 206.265)


def test_pixel_scale_falls_back_to_the_profile_optics_and_applies_binning():
    profile = get_profile("cassa8")
    unbinned = profile.get_pixel_scale(header())
    assert unbinned == pytest.approx(0.598, abs=0.001)
    assert profile.get_pixel_scale(header(XBINNING=2)) == pytest.approx(2 * unbinned)


def test_pixel_scale_from_an_existing_wcs_accounts_for_rotation():
    """CD1_1 alone is scale*cos(theta); a rotated field must not read low."""
    scale_deg, angle = 0.598 / 3600.0, np.radians(30.0)
    frame = header(
        CTYPE1="RA---TAN", CTYPE2="DEC--TAN", CRVAL1=10.0, CRVAL2=20.0,
        CRPIX1=512, CRPIX2=512,
        CD1_1=-scale_deg * np.cos(angle), CD1_2=scale_deg * np.sin(angle),
        CD2_1=scale_deg * np.sin(angle), CD2_2=scale_deg * np.cos(angle),
    )
    assert get_profile("generic").get_pixel_scale(frame) == pytest.approx(0.598, abs=0.002)


def test_pixel_scale_admits_when_nothing_can_say():
    assert get_profile("generic").get_pixel_scale(header()) is None


def test_linearity_is_the_identity_until_a_profile_measures_it():
    """The honest default: record that no correction was applied, do not invent one."""
    data = np.array([1.0, 2.0, 3.0])
    assert np.array_equal(get_profile("generic").apply_linearity(data, header()), data)


def test_a_real_frame_now_resolves_a_scale_that_phase_two_can_use():
    """The workshop frames carry SECPIX; PIXSCALE used to come out 0.0 anyway."""
    import glob
    frames = glob.glob("workshop/raw/*.fits")
    if not frames:
        pytest.skip("workshop dataset not present")
    frame_header = fits.getheader(sorted(frames)[0])
    assert get_profile("generic").get_pixel_scale(frame_header) == pytest.approx(0.4)
