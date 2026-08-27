import pytest

from cassa_photometry.config import load_config


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
