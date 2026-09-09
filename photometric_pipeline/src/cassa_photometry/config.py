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
class StepToggles:
    """Which operations a phase performs.

    Turning a step off is a legitimate thing to want -- reducing frames that
    were flat-fielded elsewhere, skipping cosmic-ray rejection on a short
    exposure, producing a catalog without calibrating it. Doing it from config
    means a user does not have to fork the code and then re-merge it forever.

    Every phase records which steps it skipped, in its log and in the product
    header, so a partially-reduced file says so rather than looking finished.
    """

    #: Steps to switch off, as a list. Equivalent to setting each to ``false``;
    #: it exists because naming five steps on one line reads better than five
    #: separate booleans, and because it pairs naturally with ``order``.
    exclude: list = field(default_factory=list)
    #: Explicit execution order. Empty means the phase's canonical order. Must
    #: list exactly the steps that run; an order that violates a declared
    #: dependency is refused at load time. See ``cassa_photometry.steps``.
    order: list = field(default_factory=list)
    #: Your own steps, as ``name: path/to/file.py:callable``.
    custom: dict = field(default_factory=dict)

    #: Field names that configure the plan rather than naming a step.
    _CONTROL = ("exclude", "order", "custom")

    @classmethod
    def step_names(cls):
        """The steps this phase declares, in declaration order."""
        return [f.name for f in fields(cls) if f.name not in cls._CONTROL]

    def enabled(self, name):
        """Whether a step runs. Unknown names are enabled, not silently off.

        ``exclude`` is honoured here as well as the booleans, so every consumer
        gets both forms without having to know about either.
        """
        if name in (self.exclude or []):
            return False
        return bool(getattr(self, name, True))

    def skipped(self):
        """The names of the steps that are switched off, by either form."""
        return sorted(
            name for name in self.step_names() if not self.enabled(name)
        )


@dataclass
class Phase1Steps(StepToggles):
    """Steps in instrument signature removal."""

    overscan: bool = True
    bias: bool = True
    dark: bool = True
    flat: bool = True
    linearity: bool = True
    cosmic_rays: bool = True
    bad_pixel_mask: bool = True
    measure_fwhm: bool = True


@dataclass
class Phase2Steps(StepToggles):
    """Steps in integration."""

    align: bool = True
    subtract_background: bool = True
    stack: bool = True
    solve_wcs: bool = True
    measure_fwhm: bool = True
    visual_qa: bool = True


@dataclass
class Phase3Steps(StepToggles):
    """Steps in photometry."""

    zero_point: bool = True
    aperture_correction: bool = True
    flux_calibration: bool = True
    catalog: bool = True
    psf_photometry: bool = True
    classification: bool = True


@dataclass
class Phase4Steps(StepToggles):
    """Stages in the diagnostics report.

    Named for the products they inspect, matching the ``stage_*`` keys in the
    report, ``metrics.json`` and the PNG filenames -- so a user reading a report
    already knows what to switch off.
    """

    stage_0_raw: bool = True
    stage_1_calibrated: bool = True
    stage_2_master: bool = True
    stage_3_photometry: bool = True


@dataclass
class Phase1Config:
    """Instrument Signature Removal parameters."""

    #: Which steps run. See ``docs/CUSTOMIZING.md``.
    steps: Phase1Steps = field(default_factory=Phase1Steps)

    # Sigma-clipping threshold when combining calibration frames. A median is
    # robust to a few outliers but not immune, and nothing else in phase 1 looks
    # at calibration frames for cosmic rays or satellite trails. 0 disables it.
    master_sigma_clip: float = 3.0
    # Detector non-linearity correction, applied in raw ADU before any other
    # step. Only does anything when the instrument profile supplies a curve --
    # the response is 1-5% non-linear near full well on a typical CMOS sensor,
    # and because the error grows with signal it tilts the magnitude scale
    # rather than offsetting it.
    apply_linearity: bool = True
    # How far a calibration frame's sensor temperature may differ from the
    # science frames' before it is reported. Dark current roughly doubles every
    # 6-7 C, so a mismatched dark removes the wrong thermal signal.
    calibration_temp_tolerance_c: float = 3.0
    # Cosmic-ray rejection (astroscrappy.detect_cosmics)
    cr_sigclip: float = 4.5
    cr_sigfrac: float = 0.3
    cr_objlim: float = 5.0
    # Largest share of a frame the cosmic-ray mask may claim before it is thrown
    # away as implausible. Cosmic rays arrive at a rate set by the sky rather
    # than by the field -- of order a few events per cm^2 per minute, marking
    # well under 0.1% of the pixels even on a large sensor in a long exposure --
    # so a mask covering percent of the frame is not cosmic rays. In practice it
    # means the fine-structure model has met a crowded field or a large resolved
    # object, whose stellar peaks are exactly the sharp features it hunts.
    # Set to 1.0 to accept any mask.
    cr_max_fraction: float = 0.01
    # Saturation threshold in raw ADU used to flag the DQ SATURATED bit. Only a
    # fallback: a profile that knows its detector's full well (from the header or
    # its hardware table) overrides this via InstrumentProfile.get_saturation.
    saturation_adu: float = 50000.0
    # Last-resort detector constants, used only when neither the header nor the
    # instrument profile can supply them. Reaching these means the error budget
    # is built on an assumption, so phase 1 logs a warning naming the frame.
    fallback_gain: float = 1.0
    fallback_read_noise: float = 10.0
    # Reduce frames that report having been calibrated already (a non-empty
    # CALSTAT, or data in electrons). Off by default: calibrating a reduced
    # frame a second time is silently wrong rather than loud.
    allow_precalibrated: bool = False
    # --- Bad-pixel mask -------------------------------------------------------
    # Each test is relative to its own master's median, so the defaults carry
    # across detectors without retuning. Set a factor to 0 to disable that test.
    # Sensitivity defects, as fractions of the flat's own LOCAL level -- the
    # flat divided by a median-smoothed copy of itself. Measuring against the
    # global median instead would condemn the corners of any strongly vignetted
    # optical train, since vignetting is optics rather than a detector defect.
    bpm_flat_low: float = 0.5
    bpm_flat_high: float = 1.5
    # Box size of that smoothing, in pixels. 0 disables it and restores the
    # older behaviour of thresholding against the normalised flat directly.
    bpm_flat_smooth_px: int = 51
    # Hot pixels: dark current above this multiple of the median dark rate...
    bpm_dark_rate_factor: float = 20.0
    # ...or this many robust sigmas above it, whichever cut is higher. The sigma
    # term is what stops a near-zero median from flagging the whole frame.
    bpm_dark_sigma: float = 5.0
    # Unstable pixels: bias frame-to-frame scatter above this multiple of the
    # median scatter, or this many robust sigmas above it, whichever is higher.
    bpm_bias_noise_factor: float = 5.0
    bpm_bias_noise_sigma: float = 5.0


@dataclass
class Phase2Config:
    """Alignment / stacking / astrometry parameters."""

    steps: Phase2Steps = field(default_factory=Phase2Steps)

    align_detection_sigma: float = 1.5
    align_min_area: int = 4
    scale_min: float = 0.65
    scale_max: float = 1.5
    background_box: tuple = (50, 50)
    background_filter: tuple = (3, 3)
    stack_sigma: float = 3.0
    stack_maxiters: int = 3
    # Frames >= this count are sigma-clipped when combining. Below it nothing is
    # rejected: with a handful of frames there is not enough information to
    # identify an outlier without discarding most of the signal.
    sigma_clip_min_frames: int = 5
    # Fewest aligned frames worth stacking at all.
    min_frames_to_stack: int = 3
    # Exclude pixels phase 1 flagged (saturated / bad / cosmic ray) from the
    # combine instead of averaging them in as good data.
    mask_bad_pixels: bool = True
    # Stacking weights. "point_source" uses 1/(sigma*FWHM)^2, optimal for stars;
    # "extended" uses 1/sigma^2, optimal for extended flux.
    stack_weight: str = "point_source"
    # Anchor ranking: SNR / FWHM**power. 0 restores ranking on SNR alone.
    anchor_sharpness_power: float = 1.0
    # Reject frames whose FWHM exceeds this multiple of the group median. 0 off.
    fwhm_reject_factor: float = 1.6
    # Aperture used to measure the frame-to-frame flux scale, in units of the
    # WORSE of the two frames' FWHM, with a fallback when neither reports one.
    # Measured on a synthetic 1.5x seeing change, where the truth is 1.000:
    #   peak pixel (the old method)  2.250   <- reads seeing as transparency
    #   aperture at 1.5x FWHM        1.097
    #   aperture at 2.0x FWHM        1.048
    #   aperture at 3.0x FWHM        1.012
    # A Moffat's wings are why this needs to be generous.
    scale_aperture_factor: float = 3.0
    scale_aperture_fwhm_default: float = 4.5
    # How frames are binned in time. "night" keeps separate observing nights
    # apart, binned on LOCAL noon so a night crossing UT midnight stays whole.
    # "none" reproduces the historical grouping exactly (all dates together).
    # A duration such as "6h" gives sub-night bins -- needed for fast variables,
    # for which nightly binning averages away the variability being measured.
    epoch_bin: str = "night"
    # Subtract a 2D background from each frame before stacking. On by default
    # (it is what makes registration and stacking well behaved), but it removes
    # real extended emission with the sky, so extended-source work may want it
    # off. The removed level is recorded as BKGLEVEL either way.
    subtract_background: bool = True
    warp_order: int = 3
    # --- Plate solving --------------------------------------------------------
    # Which backend solves the WCS. "auto" tries them in the order ASTAP,
    # solve-field, astrometry-py and uses the first that can actually run;
    # naming one forces it. ASTAP is preferred because it is a sub-megabyte
    # binary with no Python dependency, packaged for every platform this
    # pipeline supports -- including ARM Linux, where neither of the other two
    # exists -- and because it solved the CASSA test frames in 0.1 s against
    # 41 s for the in-process solver.
    solver: str = "auto"
    # Path to the ASTAP executable. None -> look for astap_cli then astap on PATH.
    astap_path: str | None = None
    # Where ASTAP star tiles live. None -> CASSA_ASTAP_DB, then
    # ~/.cache/cassa-photometry/astap. Mirrors the astrometry.net index cache.
    astap_db_dir: str | None = None
    # Which ASTAP series to fetch from. "auto" picks one from the frame's own
    # field height; naming a series forces it. d50 covers 6 deg down to well
    # under its documented 0.2 deg floor (verified solving a 0.17 deg field),
    # which spans every CASSA telescope, and is the only series published as a
    # ZIP -- the format that makes per-tile fetching possible at all.
    astap_db_series: str = "auto"
    # Where the published archives live. Only series distributed as a ZIP can be
    # fetched incrementally, which is why this points at the star_databases tree.
    astap_db_url: str = ("https://sourceforge.net/projects/astap-program/files/"
                         "star_databases/")
    # Fetch missing tiles on demand. False selects but never downloads, and
    # reports which tiles a field needs -- the air-gapped case.
    astap_db_download: bool = True
    # How many tiles either side of the field centre to keep. 1 covers a field
    # sitting on a tile edge and a pointing that is a few arcminutes out.
    astap_tile_neighbours: int = 1
    # Astrometry index directory holding a full local set. None -> resolve from
    # the CASSA_ASTROMETRY_INDEX env var, then ./astrometry_data. A populated
    # directory still wins over the on-demand cache, so existing installations
    # keep working unchanged.
    astrometry_index_dir: str | None = None

    # --- On-demand astrometry indexes -----------------------------------------
    # A field needs a handful of index files, not the whole 5 GB set. These
    # control selecting and fetching only what a frame actually requires.
    #
    # Base URL of the index series. The public Astrometry.net data server is the
    # default, so a fresh clone works with no configuration and the observatory
    # hosts nothing.
    index_url: str = "https://data.astrometry.net/"
    # Optional bearer token, for a private mirror. Nothing needs it by default.
    index_token: str | None = None
    # None -> ~/.cache/cassa-photometry/astrometry
    index_cache_dir: str | None = None
    index_cache_gb: float = 20.0
    # Search cone around the pointing when selecting index tiles, and the wider
    # one used by the fallback pass after a first solve fails.
    index_search_radius_deg: float = 3.0
    index_blind_radius_deg: float = 15.0
    # Quad-scale window to accept, as a fraction of the field size.
    # Astrometry.net's guidance is 10%-100%, but the finest indexes are by far
    # the largest files (index-4200-14 alone is 410 MB against 55 MB for
    # index-4203-14) and are rarely what solves a frame. So the first pass asks
    # for 30%-100%, which is 246 MB rather than 873 MB for a CASSA field, and
    # the widened fallback pass drops to the full 10% only if that fails.
    index_scale_lo_frac: float = 0.30
    index_scale_lo_frac_wide: float = 0.10
    index_scale_hi_frac: float = 1.00
    # Download only when this is true; otherwise selection still runs and a
    # missing file is reported with the cassa-index-fetch command that gets it.
    index_download: bool = True

    # --- Solver hints ---------------------------------------------------------
    # Half-width of the pixel-scale window handed to the solver, as a fraction.
    # Deliberately generous: a header's stated scale is usually right but can be
    # badly wrong (the iTelescope frames claim SECPIX=0.4 where the truth is
    # 0.591), and a wrong scale makes a solve FAIL where no scale only makes it
    # slow. The fallback pass widens this further before giving up on it.
    solve_scale_tolerance: float = 0.25
    solve_scale_tolerance_wide: float = 1.0
    # Warn when the solved plate scale disagrees with the header's claim by more
    # than this fraction. A permanent data-quality check: it is how a wrong
    # SECPIX gets noticed instead of silently breaking every future solve.
    solve_scale_warn_frac: float = 0.05
    solve_timeout_s: float = 120.0
    solve_timeout_wide_s: float = 600.0
    # Where the raw frames live, for the QA PDF's "true raw" panel only. Normally
    # unset: phase 1 stamps RAWDIR into every calibrated header, so the raw tree
    # is found even when it sits outside the work directory. Set this when the
    # raws have since moved, or for frames calibrated before RAWDIR existed.
    raw_dir: str | None = None


@dataclass
class Phase3Config:
    """Photometry / zero-point parameters."""

    steps: Phase3Steps = field(default_factory=Phase3Steps)

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
    # The background annulus must sit OUTSIDE the PSF wings. At 3-4x FWHM a
    # Moffat still puts ~1.25% of the star's own light in the annulus, which is
    # then subtracted off as if it were sky.
    annulus_in_factor: float = 5.0
    annulus_out_factor: float = 8.0
    zp_sigma_clip: float = 3.0
    detect_npixels: int = 5
    deblend_nlevels: int = 32
    deblend_contrast: float = 0.001
    # Objects rounder than this are treated as stars (verification tool).
    ellipticity_star_max: float = 0.15
    # PSF-fitted photometry, the matched filter for a point source. Also what
    # star/galaxy classification is based on, via MAG_PSF - MAG_AUTO.
    psf_photometry: bool = True
    # Persist every successful reference-catalog query, so a field reduced once
    # can be re-reduced with no network. This is what makes the pipeline usable
    # away from the university.
    catalog_cache: bool = True
    catalog_cache_dir: str | None = None
    # Refresh a cached query older than this many days. 0 keeps it forever.
    catalog_cache_days: float = 180.0
    # Never touch the network: use the cache only, and fail clearly when a field
    # is not in it rather than hanging on a series of timeouts.
    offline: bool = False
    # Zero-flux of the AB system, in Jy.
    ab_flux_zero_jy: float = 3631.0
    # Apply the band's AB-minus-Vega offset in the Jy conversion. B/V/R/I are
    # calibrated against Vega-based magnitudes but the conversion is the AB
    # relation, so without this B fluxes are ~8% wrong.
    apply_ab_offset: bool = True
    # Measure an aperture correction from a curve of growth, so the zero point
    # is a TOTAL-flux zero point rather than an aperture one. Without it the ZP
    # and the catalog magnitudes it calibrates sit on different flux scales, and
    # the error is a brightness-dependent tilt rather than an offset.
    aperture_correction: bool = True
    # Radius treated as "total" in the curve of growth, in units of the FWHM.
    apcor_total_factor: float = 5.0
    # Reject calibrator stars with a saturated/bad/cosmic-ray pixel inside the
    # aperture. The brightest catalog stars saturate first and are exactly the
    # ones inverse-variance weighting trusts most.
    zp_reject_flagged: bool = True
    # Keep sources with non-positive flux in the catalog, reported with their
    # uncertainty and a limiting magnitude, instead of dropping them. Dropping
    # them biases faint number counts and makes upper limits impossible.
    keep_negative_flux: bool = True


@dataclass
class Phase4Config:
    """Diagnostics parameters.

    Phase 4 measures rather than decides, so these only steer what it looks at.
    The PSF settings are shared with the FWHM measurement the science phases now
    depend on, which is why they live in config rather than as call-site
    defaults.
    """

    steps: Phase4Steps = field(default_factory=Phase4Steps)

    # PSF / FWHM estimation (phase4_diagnostics.psf.estimate_fwhm).
    fwhm_guess: float = 3.5
    detection_threshold: float = 5.0
    max_stars: int = 25
    cutout: int = 15
    # A pixel within this fraction of the frame maximum counts as near-saturated
    # in the stage-0/1 saturation map. A crude proxy, deliberately independent of
    # the DQ SATURATED bit so the two can be compared.
    near_saturation_frac: float = 0.95
    # Depth reporting.
    limiting_snr: float = 5.0
    mag_binsize: float = 0.5
    hist_bins: int = 30


@dataclass
class DetectorOverride:
    """Describe an observing setup without writing an instrument profile.

    This is the answer for the common case: a camera that does not record its
    own gain and read noise, on a telescope the pipeline ships no profile for.
    Set what you know here and reduce with ``instrument: generic``; writing a
    profile class (see ``examples/profiles/``) is only needed when the answer
    depends on the header -- a CMOS gain curve, a multi-amplifier readout.

    Values fill in **where the header is silent**, preserving the pipeline's
    order of authority: header first, then this, then the profile, then the
    configured fallbacks. Set ``override_header`` to invert the first two, which
    is occasionally necessary when acquisition software writes a wrong card --
    and is logged loudly, because a header that disagrees with the config is
    worth knowing about.
    """

    gain: float | None = None                 # e-/ADU
    read_noise: float | None = None           # e-
    saturation_adu: float | None = None       # raw ADU at which pixels clip
    pixel_scale_arcsec: float | None = None   # unbinned arcsec/pixel
    #: Extra raw-filter -> stacking-label entries, merged over the profile's map.
    filter_map: dict = field(default_factory=dict)
    #: Extra raw-filter -> science-band entries (R/G/B/V/I), merged likewise.
    science_bands: dict = field(default_factory=dict)
    #: Filters to exclude from photometric calibration entirely.
    uncalibrated_filters: list = field(default_factory=list)
    #: filter -> stand-in filter, for a setup that reuses one flat for another.
    flat_proxies: dict = field(default_factory=dict)
    #: Let these values beat the frame header rather than fill in behind it.
    override_header: bool = False

    def is_empty(self):
        """True when nothing was configured, so the profile is used unchanged."""
        return all(
            getattr(self, f.name) in (None, {}, [], False)
            for f in fields(self)
        )


@dataclass
class PipelineConfig:
    """Top-level configuration aggregating the per-phase settings."""

    phase1: Phase1Config = field(default_factory=Phase1Config)
    phase2: Phase2Config = field(default_factory=Phase2Config)
    phase3: Phase3Config = field(default_factory=Phase3Config)
    phase4: Phase4Config = field(default_factory=Phase4Config)
    #: Detector constants for a setup with no shipped profile; see
    #: :class:`DetectorOverride` and ``docs/CUSTOMIZING.md``.
    detector: DetectorOverride = field(default_factory=DetectorOverride)
    #: Instrument-profile name; see ``cassa_photometry.instruments.registry``.
    #: A ``--instrument`` flag on any command overrides this.
    instrument: str = "generic"
    #: Optional ``targets.yaml`` describing what is being observed. Frames that
    #: carry TARGNAME/OBJTYPE need no file; this retrofits archival data and
    #: holds what the telescope cannot know (comparison stars, a refined
    #: position after discovery).
    targets_file: str | None = None
    #: A profile class in a local file, as ``path/to/my_profile.py:MyProfile``.
    #: Takes precedence over ``instrument``. This is how a user runs their own
    #: profile without installing anything or editing this package.
    instrument_module: str | None = None
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

    def resolve_index_cache_dir(self):
        """Where on-demand index files are cached.

        Config value > ``CASSA_INDEX_CACHE`` > ``~/.cache/cassa-photometry/astrometry``.
        Separate from :meth:`resolve_astrometry_index_dir`, which is a full local
        set the user manages themselves and which always takes precedence.
        """
        if self.phase2.index_cache_dir:
            return os.path.abspath(os.path.expanduser(self.phase2.index_cache_dir))
        env = os.environ.get("CASSA_INDEX_CACHE")
        if env:
            return os.path.abspath(os.path.expanduser(env))
        base = os.environ.get("XDG_CACHE_HOME") or os.path.join(
            os.path.expanduser("~"), ".cache"
        )
        return os.path.join(base, "cassa-photometry", "astrometry")

    def resolve_catalog_cache_dir(self):
        """Where reference-catalog queries are cached.

        Config value > ``CASSA_CATALOG_CACHE`` > ``~/.cache/cassa-photometry/catalogs``.
        """
        if self.phase3.catalog_cache_dir:
            return os.path.abspath(os.path.expanduser(self.phase3.catalog_cache_dir))
        env = os.environ.get("CASSA_CATALOG_CACHE")
        if env:
            return os.path.abspath(os.path.expanduser(env))
        base = os.environ.get("XDG_CACHE_HOME") or os.path.join(
            os.path.expanduser("~"), ".cache"
        )
        return os.path.join(base, "cassa-photometry", "catalogs")

    def resolve_index_token(self):
        """Bearer token for a private index mirror, or None."""
        return self.phase2.index_token or os.environ.get("CASSA_INDEX_TOKEN") or None


class ConfigError(ValueError):
    """A config file asked for something the pipeline cannot honour."""


#: Type names understood in an annotation for a field that defaults to None.
_OPTIONAL_TYPES = (("float", float), ("int", int), ("str", str))


def _coerce_optional(path, annotation, value):
    """Validate a field whose default is ``None``, using its annotation.

    A ``None`` default carries no type information, so these fields -- the
    optional detector constants and paths -- used to accept anything at all,
    which made ``detector: {gain: "not-a-number"}`` a silent no-op. The
    annotation is a *string* here (this module uses postponed annotations), so
    it is matched textually rather than resolved.
    """
    annotation = str(annotation or "")
    if value is None:
        return None
    for name, expected in _OPTIONAL_TYPES:
        if name not in annotation:
            continue
        if expected is float:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ConfigError(f"{path}: expected a number, got {value!r}")
            return float(value)
        if expected is int:
            if isinstance(value, bool) or not isinstance(value, int):
                raise ConfigError(f"{path}: expected an integer, got {value!r}")
            return int(value)
        if not isinstance(value, str):
            raise ConfigError(f"{path}: expected a string, got {value!r}")
        return value
    return value


def _coerce(path, current, value, annotation=None):
    """Coerce a YAML value to the type of the field it overrides.

    The check is against the *default value's* type rather than the annotation,
    because this module uses postponed annotations (so ``field.type`` is a
    string like ``"str | None"``). That is also the more useful check: it is the
    default that documents what the field means. Fields that default to ``None``
    are the exception and fall back to the annotation.

    Without this, a typo such as ``phase3: {fwhm: "4.O"}`` was accepted silently
    and only surfaced much later as a confusing error inside photutils -- or,
    worse, not at all.
    """
    if current is None:
        return _coerce_optional(path, annotation, value)

    if isinstance(current, bool):
        if isinstance(value, bool):
            return value
        raise ConfigError(f"{path}: expected true/false, got {value!r}")

    if isinstance(current, int) and not isinstance(current, bool):
        if isinstance(value, bool):
            raise ConfigError(f"{path}: expected an integer, got {value!r}")
        if isinstance(value, int):
            return value
        if isinstance(value, float) and float(value).is_integer():
            return int(value)
        raise ConfigError(f"{path}: expected an integer, got {value!r}")

    if isinstance(current, float):
        if isinstance(value, bool):
            raise ConfigError(f"{path}: expected a number, got {value!r}")
        if isinstance(value, (int, float)):
            return float(value)
        raise ConfigError(f"{path}: expected a number, got {value!r}")

    if isinstance(current, str):
        if isinstance(value, str):
            return value
        raise ConfigError(f"{path}: expected a string, got {value!r}")

    if isinstance(current, tuple):
        if isinstance(value, (list, tuple)):
            return tuple(value)
        raise ConfigError(f"{path}: expected a list, got {value!r}")

    if isinstance(current, dict):
        if isinstance(value, dict):
            return dict(value)
        raise ConfigError(f"{path}: expected a mapping, got {value!r}")

    if isinstance(current, list):
        if isinstance(value, (list, tuple)):
            return list(value)
        raise ConfigError(f"{path}: expected a list, got {value!r}")

    return value


def _merge(dc, overrides, prefix=""):
    """Recursively apply a plain dict of ``overrides`` onto a dataclass instance.

    Unknown keys and type mismatches are both errors, and both name the full
    dotted path so the message points at the line in the YAML file.
    """
    if not overrides:
        return dc
    if not isinstance(overrides, dict):
        raise ConfigError(f"{prefix or 'config'}: expected a mapping, got {overrides!r}")

    valid = {f.name: f for f in fields(dc)}
    for key, value in overrides.items():
        path = f"{prefix}{key}"
        if key not in valid:
            raise KeyError(
                f"Unknown config key: {path!r}. "
                f"Valid keys here: {', '.join(sorted(valid))}."
            )
        current = getattr(dc, key)
        if is_dataclass(current):
            if not isinstance(value, dict):
                raise ConfigError(f"{path}: expected a mapping of settings, got {value!r}")
            _merge(current, value, prefix=f"{path}.")
        else:
            setattr(dc, key, _coerce(path, current, value, valid[key].type))
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

    _fold_legacy_step_flags(config)
    _validate_step_plans(config)

    # log_level was declared but never read; applying it here is what makes
    # `log_level: DEBUG` in a config file reach every phase's logger.
    from cassa_photometry.logging_utils import set_default_level

    set_default_level(config.log_level)
    return config


#: Settings that duplicate a ``steps`` toggle, as
#: (phase, legacy field, step name). Each pair used to be ANDed at the call
#: site, so switching one off left the other still blocking and the user could
#: not tell which. The ``steps`` toggle is now the single authority and these
#: fold into it.
_LEGACY_STEP_FLAGS = [
    ("phase1", "apply_linearity", "linearity"),
    ("phase2", "subtract_background", "subtract_background"),
    ("phase3", "psf_photometry", "psf_photometry"),
    ("phase3", "aperture_correction", "aperture_correction"),
]


def _validate_step_plans(config):
    """Resolve every phase's plan now, so a bad one fails before any work.

    A misordered plan discovered three hours into a reduction is a wasted night.
    Resolving here also loads any custom step, so a typo in its path is reported
    against the config file rather than surfacing much later as an ImportError
    from inside a phase.
    """
    from cassa_photometry.steps import resolve_plan

    for phase in ("phase1", "phase2", "phase3", "phase4"):
        resolve_plan(phase, getattr(config, phase).steps)


def _fold_legacy_step_flags(config):
    """Let a legacy ``false`` still switch its step off, and say so.

    Kept working rather than removed: these names appear in existing config
    files, in ``README.md`` and in the workshop notebooks, and a setting that
    silently stops having an effect is the same defect this whole change exists
    to fix. Setting one to ``false`` therefore still disables the step; it is
    only the *duplication* that is gone.
    """
    from cassa_photometry.logging_utils import get_logger

    for phase, legacy, step in _LEGACY_STEP_FLAGS:
        section = getattr(config, phase)
        if getattr(section, legacy):
            continue
        get_logger("cassa_config").warning(
            "%s.%s is deprecated; use %s.steps.%s instead. Honouring it: %s "
            "is switched off.", phase, legacy, phase, step, step,
        )
        setattr(section.steps, step, False)
