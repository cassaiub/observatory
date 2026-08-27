"""Central configuration for the sensor-characterization pipeline.

Every value that used to be hardcoded across the two scripts lives here as a
dataclass field. Defaults reproduce the original QHYCCD miniCAM8M / IMX585
reference campaign; a YAML file overrides any subset.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


@dataclass
class SensorConfig:
    """Detector format / range."""

    image_size: tuple = (500, 500)
    max_adu: int = 65535
    base_bias_adu: float = 500.0


@dataclass
class CaptureConfig:
    """The capture grid the campaign sweeps over (also used by the simulator)."""

    driver_gains: list = field(default_factory=lambda: [0, 50, 100, 150, 200])
    exposure_min_s: float = 0.1
    exposure_max_s: float = 12.0
    exposure_steps: int = 20
    flat_electrons_per_s: float = 5000.0                 # simulator flat brightness
    temperatures_c: list = field(default_factory=lambda: [-20, -15, -10, -5, 0, 5, 10, 15, 20])
    dark_temp_exposure_s: float = 300.0
    linearity_temp_c: float = 0.0
    linearity_exposures_s: list = field(default_factory=lambda: [10, 30, 60, 120, 300, 600])


@dataclass
class AnalysisConfig:
    """Analysis parameters (robust statistics + PTC region + spectral grid)."""

    ptc_linear_lo: float = 0.1        # PTC linear-fit region: fraction of max mean signal
    ptc_linear_hi: float = 0.7
    roi_fraction: float = 0.5         # central ROI side fraction for frame statistics
    sigma_clip: float = 5.0           # reject hot / telegraph pixels
    wavelength_min_nm: float = 350.0
    wavelength_max_nm: float = 1000.0
    wavelength_points: int = 500


@dataclass
class PathsConfig:
    data_dir: str = "./data"
    out_dir: str = "./camchar_output"


@dataclass
class CamCharConfig:
    """Top-level configuration aggregating all groups."""

    sensor: SensorConfig = field(default_factory=SensorConfig)
    capture: CaptureConfig = field(default_factory=CaptureConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    paths: PathsConfig = field(default_factory=PathsConfig)
    log_level: str = "INFO"


def _merge(dc, overrides):
    """Recursively apply a plain dict of overrides onto a dataclass instance."""
    if not overrides:
        return dc
    valid = {f.name: f for f in fields(dc)}
    for key, value in overrides.items():
        if key not in valid:
            raise KeyError(f"Unknown config key: {key!r}")
        current = getattr(dc, key)
        if is_dataclass(current) and isinstance(value, dict):
            _merge(current, value)
        elif isinstance(current, tuple) and isinstance(value, list):
            setattr(dc, key, tuple(value))
        else:
            setattr(dc, key, value)
    return dc


def load_config(path=None):
    """Build a :class:`CamCharConfig`, optionally overriding defaults from YAML."""
    config = CamCharConfig()
    if path:
        if yaml is None:  # pragma: no cover
            raise ImportError("PyYAML is required to load a config file.")
        with open(path) as handle:
            data = yaml.safe_load(handle) or {}
        _merge(config, data)
    return config
