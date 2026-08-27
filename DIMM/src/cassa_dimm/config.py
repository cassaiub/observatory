"""Central configuration for the DIMM seeing monitor.

Every tunable (optical constants, site coordinates, watch-folder behaviour,
detection thresholds, QC limits, output paths) lives here as a dataclass field.
Defaults reproduce the 8-inch SkyWatcher / EQ6R-Pro reference setup; a YAML file
overrides any subset.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

_ARCSEC_TO_RAD = 4.84813681109536e-6


@dataclass
class HardwareConfig:
    """Optical / detector constants of the DIMM instrument."""

    pixel_size_um: float = 3.75         # camera pixel pitch
    focal_length_mm: float = 1000.0     # 8-inch f/5 Newtonian
    aperture_diameter_m: float = 0.06   # D: diameter of one mask hole
    hole_separation_m: float = 0.15     # d: centre-to-centre hole separation
    wavelength_m: float = 500e-9        # reference wavelength
    gain: float = 1.0                   # e-/ADU (for SNR estimates)
    read_noise: float = 10.0            # e-
    saturation_adu: float = 60000.0
    # Optional fixed prism doublet vector in pixels; None -> auto-detect from data.
    prism_dx_pix: "float | None" = None
    prism_dy_pix: "float | None" = None

    @property
    def plate_scale_arcsec(self):
        """Plate scale in arcsec/pixel."""
        return 206.265 * self.pixel_size_um / self.focal_length_mm

    @property
    def scale_rad(self):
        """Plate scale in radians/pixel."""
        return self.plate_scale_arcsec * _ARCSEC_TO_RAD


@dataclass
class SiteConfig:
    """Observatory location, used to compute target altitude / airmass."""

    latitude_deg: float = 0.0     # SET to the CASSA site latitude
    longitude_deg: float = 0.0    # SET to the CASSA site longitude (East positive)
    elevation_m: float = 0.0


@dataclass
class WatchConfig:
    """Watch-folder / rolling-window behaviour."""

    folder: str = "./dimm_watch"
    pattern: str = "*.fit*"
    poll_interval_s: float = 2.0       # how often to scan for new frames
    window_frames: int = 150           # rolling buffer size (frames per estimate)
    cadence_frames: int = 50           # emit an estimate every this many new frames
    stabilization_s: float = 0.5       # a file must be size-stable this long before ingest
    process_existing: bool = False     # also process files already in the folder at startup


@dataclass
class DetectionConfig:
    """Source detection and doublet pairing."""

    snr_threshold: float = 5.0
    fwhm_pixels: float = 3.5
    centroid_box: int = 9              # odd box side for windowed centroiding
    max_sources: int = 60
    pair_tol_pixels: float = 2.0       # tolerance when matching the prism vector
    edge_margin_pixels: int = 6
    ellipticity_box: int = 11          # box for the per-source ellipticity moment
    # Wide-field controls (0 / defaults reproduce the small-field behaviour):
    central_radius_arcmin: float = 0.0  # >0: use only spots within this radius of centre
    ema_window: int = 50               # frames averaged into the streaming reference
    warmup_frames: int = 25            # frames before the first streaming detection
    refresh_frames: int = 300          # re-detect the doublet geometry every N frames


@dataclass
class QCConfig:
    """Quality-control limits for accepting a window / doublet."""

    min_frames: int = 30
    min_valid_fraction: float = 0.5    # fraction of frames a doublet must be measured in
    min_doublets: int = 1
    sigma_clip: float = 4.0            # frame-level outlier rejection on the differential motion
    max_lt_discrepancy: float = 0.6    # max |eps_l - eps_t| / mean before flagging
    max_ellipticity: float = 0.4       # spot elongation (wind-shake indicator)


@dataclass
class OutputConfig:
    """Where results are written."""

    log_dir: str = "./dimm_output"
    csv_name: str = "seeing_log.csv"
    jsonl_name: str = "seeing_log.jsonl"
    status_name: str = "status_latest.json"
    plot_name: str = "seeing_timeseries.png"
    plot_max_points: int = 720         # points kept on the rolling live plot


@dataclass
class DimmConfig:
    """Top-level configuration aggregating all groups."""

    hardware: HardwareConfig = field(default_factory=HardwareConfig)
    site: SiteConfig = field(default_factory=SiteConfig)
    watch: WatchConfig = field(default_factory=WatchConfig)
    detection: DetectionConfig = field(default_factory=DetectionConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
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
        else:
            setattr(dc, key, value)
    return dc


def load_config(path=None):
    """Build a :class:`DimmConfig`, optionally overriding defaults from a YAML file."""
    config = DimmConfig()
    if path:
        if yaml is None:  # pragma: no cover
            raise ImportError("PyYAML is required to load a config file.")
        with open(path) as handle:
            data = yaml.safe_load(handle) or {}
        _merge(config, data)
    return config
