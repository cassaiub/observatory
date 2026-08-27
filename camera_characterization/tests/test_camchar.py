import os
import numpy as np
import pytest

from cassa_camchar.config import load_config
from cassa_camchar.io import robust_diff_variance
from cassa_camchar.simulate import simulate_campaign
from cassa_camchar.pipeline import analyze


def _fast_config():
    cfg = load_config()
    cfg.sensor.image_size = (200, 200)
    cfg.capture.driver_gains = [0, 100]
    cfg.capture.exposure_steps = 10
    cfg.capture.temperatures_c = [-20, -10, 0, 10, 20]
    cfg.analysis.wavelength_points = 100
    return cfg


def test_sim_analyze_roundtrip(tmp_path):
    """Simulate a campaign with a known model and recover it."""
    cfg = _fast_config()
    data, out = tmp_path / "data", tmp_path / "out"
    simulate_campaign(str(data), config=cfg, seed=1)
    res = analyze(str(data), config=cfg, out_dir=str(out))

    i0 = res.driver_gains.index(0)
    assert abs(res.system_gain[i0] - 4.0) < 0.4        # true gain 4.0 e-/ADU at driver 0
    assert abs(res.read_noise_e[i0] - 4.5) < 1.0       # true read noise 4.5 e-
    assert res.dark_current[-1] > res.dark_current[0]  # rises with temperature

    assert os.path.exists(out / "characterization_results.csv")
    assert os.path.exists(out / "characterization_results.json")
    assert os.path.exists(out / "CASSA_Sensor_Characterization_Report.png")


def test_roi_sigma_clip_rejects_hot_pixel():
    rng = np.random.default_rng(0)
    a = rng.normal(1000, 10, (200, 200))
    b = rng.normal(1000, 10, (200, 200))
    v_clean = robust_diff_variance(a, b, 0.5, 5.0)

    a_hot = a.copy()
    a_hot[100, 100] = 60000.0                          # one hot pixel
    v_clipped = robust_diff_variance(a_hot, b, 0.5, 5.0)
    v_noclip = robust_diff_variance(a_hot, b, 0.5, 1e9)  # effectively no clipping

    assert abs(v_clipped - v_clean) / v_clean < 0.1    # sigma-clip removes it
    assert v_noclip > 2 * v_clean                      # unclipped, it dominates


def test_config_override(tmp_path):
    yaml_file = tmp_path / "c.yaml"
    yaml_file.write_text("analysis:\n  sigma_clip: 3.0\nsensor:\n  image_size: [256, 256]\n")
    cfg = load_config(str(yaml_file))
    assert cfg.analysis.sigma_clip == 3.0
    assert cfg.sensor.image_size == (256, 256)
    assert cfg.analysis.roi_fraction == 0.5            # default preserved


def test_spectral_synthetic_fallback(tmp_path):
    from cassa_camchar.spectral import analyze_spectra
    cfg = load_config()
    cfg.analysis.wavelength_points = 100
    spectra = analyze_spectra(str(tmp_path), cfg)       # empty dir -> synthetic
    assert spectra["qe_synthetic"] is True
    assert 85 < spectra["qe_peak_pct"] <= 92
    assert "H-Alpha" in spectra["filter_metrics"]


@pytest.mark.parametrize("func", ["analyze", "simulate"])
def test_cli_help(func, monkeypatch):
    from cassa_camchar import cli
    prog = "cassa-camchar-analyze" if func == "analyze" else "cassa-camchar-sim"
    monkeypatch.setattr("sys.argv", [prog, "--help"])
    with pytest.raises(SystemExit) as exc:
        getattr(cli, func)()
    assert exc.value.code == 0
