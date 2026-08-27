import os
import numpy as np

from cassa_dimm.config import load_config
from cassa_dimm.seeing import (
    response_coefficients, variance_to_r0, r0_to_seeing_arcsec,
    kasten_young_airmass, correct_seeing_to_zenith,
)
from cassa_dimm.detect import Source, estimate_prism_vector, pair_doublets
from cassa_dimm.simulate import simulate_frames
from cassa_dimm.io import read_cube, extract_meta
from cassa_dimm.window import estimate_window


def test_physics_roundtrip():
    """Inject the variance the reducer expects for a known r0 and recover it."""
    D, d, lam = 0.06, 0.15, 500e-9
    r0_true = 0.10
    c_l, c_t = response_coefficients(D, d)
    var_l = 2 * lam ** 2 * r0_true ** (-5 / 3) * c_l
    r0 = variance_to_r0(var_l, c_l, lam)
    assert np.isclose(r0, r0_true, rtol=1e-6)
    # seeing FWHM = 0.98 lambda / r0
    seeing = r0_to_seeing_arcsec(r0, lam)
    assert np.isclose(seeing, 0.98 * lam / r0 / 4.84813681109536e-6, rtol=1e-6)


def test_airmass():
    assert np.isclose(kasten_young_airmass(90.0), 1.0, atol=1e-3)
    assert kasten_young_airmass(30.0) > 1.9        # ~2 airmasses at 30 deg alt
    # At airmass 1, zenith seeing == observed.
    assert np.isclose(correct_seeing_to_zenith(1.5, 1.0), 1.5, rtol=1e-6)
    # Higher airmass -> smaller zenith-corrected value.
    assert correct_seeing_to_zenith(1.5, 2.0) < 1.5


def test_prism_vector_and_pairing():
    cfg = load_config()
    vec = np.array([40.0, 0.0])
    rng = np.random.default_rng(0)
    sources = []
    for _ in range(3):
        bx, by = rng.uniform(30, 160), rng.uniform(30, 160)
        sources.append(Source(bx, by, 1000, 1000, 50))
        sources.append(Source(bx + vec[0], by + vec[1], 1000, 1000, 50))
    rng.shuffle(sources)
    est, _ = estimate_prism_vector(sources, cfg)
    assert np.allclose(np.abs(est), vec, atol=1.0)
    doublets = pair_doublets(sources, est, cfg)
    assert len(doublets) == 3


def test_simulator_recovers_seeing(tmp_path):
    cfg = load_config()
    injected = 1.5
    cube_path = simulate_frames(str(tmp_path), seeing_arcsec=injected, n_stars=3,
                                n_frames=400, as_cube=True, common_motion_pix=1.0,
                                seed=1, config=cfg)
    cube, header = read_cube(cube_path)
    frames = list(cube)
    metas = [extract_meta(header)] * len(frames)
    result = estimate_window(frames, metas, cfg)
    assert result.n_stars >= 2
    assert abs(result.seeing_raw - injected) < 0.3      # recovered within ~0.3"
    # altitude=90 -> airmass ~1 -> zenith ~ raw
    assert np.isclose(result.seeing_zenith, result.seeing_raw, rtol=1e-2)


def test_stream_recovers_seeing(tmp_path):
    """Streaming path (stamp buffering, warmup, auto-detect) recovers the seeing."""
    from cassa_dimm.stream import StreamEstimator
    cfg = load_config()
    cube_path = simulate_frames(str(tmp_path), seeing_arcsec=1.5, n_stars=3,
                                n_frames=400, as_cube=True, common_motion_pix=1.0,
                                seed=3, config=cfg)
    cube, header = read_cube(cube_path)
    meta = extract_meta(header)
    est = StreamEstimator(cfg, max_window=400)
    for fr in cube:
        est.add_frame(fr, meta)                # allow_refresh=True -> self-detects
    # Stamp buffering: no full-frame window is retained, only one EMA reference.
    assert est.ref is not None and est.ref.ndim == 2
    result = est.estimate()
    assert result.n_stars >= 2
    assert abs(result.seeing_raw - 1.5) < 0.35


def test_ellipticity_measure():
    from cassa_dimm.detect import _stamp_ellipticity
    yy, xx = np.mgrid[0:21, 0:21]
    round_img = np.exp(-((xx - 10) ** 2 + (yy - 10) ** 2) / (2 * 2.0 ** 2))
    elong_img = np.exp(-((xx - 10) ** 2 / (2 * 2.0 ** 2) + (yy - 10) ** 2 / (2 * 5.0 ** 2)))
    assert _stamp_ellipticity(round_img, 10, 10, 15) < 0.1
    assert _stamp_ellipticity(elong_img, 10, 10, 15) > 0.4


def test_central_field_filter():
    from cassa_dimm.detect import detect_sources
    cfg = load_config()
    cfg.detection.central_radius_arcmin = 1.0    # ~78 px radius at the default plate scale
    cfg.detection.snr_threshold = 3.0
    rng = np.random.default_rng(0)
    ny = nx = 200
    yy, xx = np.mgrid[0:ny, 0:nx]
    img = rng.normal(100.0, 3.0, (ny, nx))
    for (sx, sy) in [(100, 100), (180, 180)]:      # centre + far corner
        img += 5000.0 * np.exp(-((xx - sx) ** 2 + (yy - sy) ** 2) / (2 * 2.0 ** 2))
    srcs = detect_sources(img, cfg)
    assert any(abs(s.x - 100) < 3 and abs(s.y - 100) < 3 for s in srcs)   # centre kept
    assert not any(abs(s.x - 180) < 3 for s in srcs)                       # corner rejected


def test_config_override(tmp_path):
    yaml_file = tmp_path / "c.yaml"
    yaml_file.write_text("hardware:\n  aperture_diameter_m: 0.08\nqc:\n  min_frames: 10\n")
    cfg = load_config(str(yaml_file))
    assert cfg.hardware.aperture_diameter_m == 0.08
    assert cfg.qc.min_frames == 10
    assert cfg.hardware.hole_separation_m == 0.15       # default preserved


def test_monitor_emit_writes_outputs(tmp_path):
    from cassa_dimm.monitor import DimmMonitor
    cfg = load_config()
    cfg.output.log_dir = str(tmp_path / "out")
    cube_path = simulate_frames(str(tmp_path / "sim"), seeing_arcsec=1.2, n_stars=2,
                                n_frames=200, as_cube=True, seed=2, config=cfg)
    cube, header = read_cube(cube_path)

    mon = DimmMonitor(config=cfg)
    meta = extract_meta(header)
    for fr in cube:
        mon.estimator.add_frame(fr, meta)
    mon._emit()

    assert os.path.exists(mon.csv_path)
    assert os.path.exists(mon.status_path)
    with open(mon.csv_path) as fh:
        lines = fh.read().strip().splitlines()
    assert len(lines) >= 2  # header + at least one record
