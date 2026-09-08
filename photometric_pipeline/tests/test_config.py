import pytest

from cassa_photometry.config import ConfigError, load_config


def test_defaults():
    cfg = load_config()
    assert cfg.phase3.fwhm == 3.5
    assert cfg.phase1.cr_sigclip == 4.5
    assert cfg.phase2.scale_min == 0.65


def test_yaml_override(tmp_path):
    yaml_file = tmp_path / "cfg.yaml"
    yaml_file.write_text("phase3:\n  fwhm: 5.0\nphase1:\n  saturation_adu: 60000\n")
    cfg = load_config(str(yaml_file))
    assert cfg.phase3.fwhm == 5.0
    assert cfg.phase1.saturation_adu == 60000
    # Untouched values keep their defaults.
    assert cfg.phase3.detection_threshold == 5.0


def test_unknown_key_raises(tmp_path):
    yaml_file = tmp_path / "bad.yaml"
    yaml_file.write_text("phase3:\n  not_a_real_key: 1\n")
    with pytest.raises(KeyError):
        load_config(str(yaml_file))


def test_astrometry_index_env(monkeypatch):
    monkeypatch.setenv("CASSA_ASTROMETRY_INDEX", "/tmp/some_index")
    cfg = load_config()
    assert cfg.resolve_astrometry_index_dir() == "/tmp/some_index"


# --- Step toggles (WP10) ------------------------------------------------------

def test_a_step_can_be_switched_off_from_config(tmp_path):
    """Turning a step off must not require forking the code."""
    path = tmp_path / "cfg.yaml"
    path.write_text("phase1:\n  steps:\n    cosmic_rays: false\n    flat: false\n")
    config = load_config(str(path))

    assert not config.phase1.steps.enabled("cosmic_rays")
    assert not config.phase1.steps.enabled("flat")
    assert config.phase1.steps.enabled("bias")
    assert config.phase1.steps.skipped() == ["cosmic_rays", "flat"]


def test_every_phase_has_step_toggles():
    config = load_config()
    for steps, expected in (
        (config.phase1.steps, "cosmic_rays"),
        (config.phase2.steps, "solve_wcs"),
        (config.phase3.steps, "zero_point"),
    ):
        assert steps.enabled(expected)
        assert steps.skipped() == []


def test_an_unknown_step_name_is_enabled_rather_than_silently_off():
    """Failing open matters: a typo must not quietly disable a reduction step."""
    assert load_config().phase1.steps.enabled("no_such_step")


def test_a_mistyped_step_value_is_rejected(tmp_path):
    path = tmp_path / "cfg.yaml"
    path.write_text("phase1:\n  steps:\n    cosmic_rays: maybe\n")
    with pytest.raises(ConfigError, match="phase1.steps.cosmic_rays"):
        load_config(str(path))


def test_an_unknown_step_key_is_rejected(tmp_path):
    path = tmp_path / "cfg.yaml"
    path.write_text("phase1:\n  steps:\n    cosmicrays: false\n")
    with pytest.raises(KeyError, match="phase1.steps.cosmicrays"):
        load_config(str(path))
