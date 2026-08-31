"""Central configuration for the pipeline.

Every value that used to be a magic number scattered across the three phases
lives here as a dataclass field. Defaults reproduce the original behaviour;
users can override any subset through a YAML file (``load_config("my.yaml")``)
or, for the astrometry index directory, the ``CASSA_ASTROMETRY_INDEX``
environment variable.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field, fields, is_dataclass

try:  # PyYAML is a hard dependency, but keep import failure legible.
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


@dataclass
class Phase1Config:
    """Instrument Signature Removal parameters."""

    # Cosmic-ray rejection (astroscrappy.detect_cosmics)
    cr_sigclip: float = 4.5
    cr_sigfrac: float = 0.3
    cr_objlim: float = 5.0
    # Saturation threshold in raw ADU used to flag the DQ SATURATED bit.
    saturation_adu: float = 50000.0
    # Bad-pixel mask thresholds, as fractions of the normalised master flat (=1.0).
    bpm_flat_low: float = 0.5
    bpm_flat_high: float = 1.5


@dataclass
class Phase2Config:
    """Alignment / stacking / astrometry parameters."""

    align_detection_sigma: float = 1.5
    align_min_area: int = 4
    scale_min: float = 0.65
    scale_max: float = 1.5
    background_box: tuple = (50, 50)
    background_filter: tuple = (3, 3)
    stack_sigma: float = 3.0
    stack_maxiters: int = 3
    # Frames >= this count use sigma-clipping; fewer use min/max rejection.
    sigma_clip_min_frames: int = 5
    warp_order: int = 3
    # Astrometry index directory. None -> resolve from env var / ./astrometry_data.
    astrometry_index_dir: str | None = None


@dataclass
class Phase3Config:
    """Photometry / zero-point parameters."""

    fwhm: float = 3.5
    detection_threshold: float = 5.0
    zp_match_tol_arcsec: float = 2.0
    # Cross-match radius for the verification tool. When the WCS carries an
    # astrometric RMS (ASTRMS, written by phase 2), the radius is derived from it
    # as sigma * ASTRMS, clamped to the min/max below; otherwise
    # ``zp_match_tol_arcsec`` above is used as a fixed fallback.
    match_radius_sigma: float = 3.0
    match_radius_min_arcsec: float = 1.0
    match_radius_max_arcsec: float = 5.0
    aperture_r_factor: float = 2.0     # aperture radius = factor * FWHM
    annulus_in_factor: float = 3.0
    annulus_out_factor: float = 4.0
    zp_sigma_clip: float = 3.0
    detect_npixels: int = 5
    deblend_nlevels: int = 32
    deblend_contrast: float = 0.001
    # Objects rounder than this are treated as stars (verification tool).
    ellipticity_star_max: float = 0.15
    # AB magnitude of a source with 1 electron/s ... actually the AB zero flux.
    ab_flux_zero_jy: float = 3631.0


@dataclass
class PipelineConfig:
    """Top-level configuration aggregating the per-phase settings."""

    phase1: Phase1Config = field(default_factory=Phase1Config)
    phase2: Phase2Config = field(default_factory=Phase2Config)
    phase3: Phase3Config = field(default_factory=Phase3Config)
    log_level: str = "INFO"

    def resolve_astrometry_index_dir(self):
        """Return the astrometry index directory using the resolution order.

        CLI/config value > ``CASSA_ASTROMETRY_INDEX`` env var > ``./astrometry_data``.
        """
        if self.phase2.astrometry_index_dir:
            return os.path.abspath(os.path.expanduser(self.phase2.astrometry_index_dir))
        env = os.environ.get("CASSA_ASTROMETRY_INDEX")
        if env:
            return os.path.abspath(os.path.expanduser(env))
        return os.path.abspath("astrometry_data")


def _merge(dc, overrides):
    """Recursively apply a plain dict of ``overrides`` onto a dataclass instance."""
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
    """Build a :class:`PipelineConfig`, optionally overriding defaults from YAML.

    Parameters
    ----------
    path : str, optional
        Path to a YAML file whose (nested) keys mirror the dataclass fields,
        e.g. ``{"phase3": {"fwhm": 4.0}}``. Unknown keys raise ``KeyError``.
    """
    config = PipelineConfig()
    if path:
        if yaml is None:  # pragma: no cover
            raise ImportError("PyYAML is required to load a config file.")
        with open(path) as handle:
            data = yaml.safe_load(handle) or {}
        _merge(config, data)
    return config
