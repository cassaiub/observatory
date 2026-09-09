"""Phase 2 orchestrator: align, stack, and astrometrically solve calibrated frames.

The stack is inverse-variance weighted and now carries a propagated ``ERR``
plane through to the master (previously the variance was computed only for a QA
log line and discarded). The arbitrary ``+100`` pedestal has been removed so the
output remains on a physical electron scale for photometry.
"""

import gc
import glob
import os

import astroalign as aa
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
from scipy.ndimage import median_filter

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import (
    DQ_BAD_PIXEL,
    DQ_COSMIC_RAY,
    DQ_FLAG_NAMES,
    DQ_NO_DATA,
    DQ_REJECTED,
    DQ_SATURATED,
)
from cassa_photometry.instruments import get_profile
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.paths import sibling_phase_dir
from cassa_photometry.phase2_integration.epochs import epoch_key
from cassa_photometry.phase2_integration.hardware import HardwareManager
from cassa_photometry.phase2_integration.io import FITSHandler
from cassa_photometry.phase2_integration.math_utils import MathEngine
from cassa_photometry.phase2_integration.models import TargetGroup
from cassa_photometry.phase2_integration.visuals import VisualQAGenerator
from cassa_photometry.phase2_integration.wcs import WCSSolver
from cassa_photometry.steps import resolve_plan


def _variance_from(err, shape, fallback_noise):
    """Return a per-pixel variance plane, falling back to a flat scalar variance."""
    if err is not None:
        return np.asarray(err, dtype=np.float64) ** 2
    return np.full(shape, float(fallback_noise) ** 2, dtype=np.float64)


def _mean_of(values):
    """Mean of the values that are actually known, or None."""
    known = [float(v) for v in values if v is not None and np.isfinite(v)]
    return float(np.mean(known)) if known else None


def _std_of(values):
    known = [float(v) for v in values if v is not None and np.isfinite(v)]
    return float(np.std(known)) if len(known) > 1 else None


def _min_of(values):
    known = [float(v) for v in values if v is not None and np.isfinite(v)]
    return float(np.min(known)) if known else None


def _max_of(values):
    known = [float(v) for v in values if v is not None and np.isfinite(v)]
    return float(np.max(known)) if known else None


def _positive_or_none(value):
    """A header value as a positive float, or None when it is absent or absurd."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) and number > 0 else None


def _as_dq(dq, shape):
    """A DQ plane as int32 of the right shape, or zeros when absent."""
    if dq is None:
        return np.zeros(shape, dtype=np.int32)
    array = np.asarray(dq)
    if array.shape != tuple(shape):
        return np.zeros(shape, dtype=np.int32)
    return array.astype(np.int32)


def _bad_pixel_mask(dq, shape):
    """Pixels that must not contribute to a combine.

    Saturated, known-bad and cosmic-ray-hit pixels carry no usable signal;
    averaging them in as good data is how a saturated stellar core ends up in a
    master looking merely bright. ``NO_DATA`` is handled by the finiteness test.
    """
    plane = _as_dq(dq, shape)
    return (plane & (DQ_SATURATED | DQ_BAD_PIXEL | DQ_COSMIC_RAY)) != 0


def _combine_dq(dq_stack, n_used, master):
    """Master DQ: the OR of the contributing frames, plus coverage.

    A bit is set when *every* contributing frame had it -- a pixel rescued by
    even one good frame is good in the stack. ``DQ_NO_DATA`` marks pixels no
    frame covered, and ``DQ_REJECTED`` those the combine threw out entirely.
    """
    cube = np.array(dq_stack, dtype=np.int32)
    n_frames = cube.shape[0]

    # A defect present in every frame survives; one present in some does not.
    always = np.bitwise_and.reduce(cube, axis=0)
    dq = always.astype(np.int32)

    no_data = (~np.isfinite(master)) | (np.asarray(n_used) <= 0)
    dq[no_data] |= DQ_NO_DATA
    partial = (np.asarray(n_used) < n_frames) & ~no_data
    dq[partial] |= DQ_REJECTED
    return dq


def _dq_summary(dq):
    """Per-flag pixel counts, for the run log."""
    return {
        name: int(np.count_nonzero(np.asarray(dq) & bit))
        for bit, name in DQ_FLAG_NAMES.items()
    }


class IntegrationPipeline:
    def __init__(self, input_dir, output_dir=None, keep_temps=False, config=None,
                 logger=None, instrument=None, assume_yes=False, cpu_level=None):
        self.input_dir = input_dir
        # Phase 2 writes to <work>/phase2 next to the phase 1 directory it reads.
        self.output_dir = output_dir or sibling_phase_dir(input_dir, 2)
        self.keep_temps = keep_temps
        self.config = config or load_config()
        self.instrument = instrument or get_profile(self.config.instrument, config=self.config)
        self.hw = HardwareManager(logger=logger, config=self.config)
        # Carried so setup() never has to decide whether prompting is possible.
        self.assume_yes = assume_yes
        self.cpu_level = cpu_level
        self.target_groups = []
        self.run_dir = None
        self.logger = logger

    def setup(self):
        cores = self.hw.allocate_resources(
            assume_yes=self.assume_yes, cpu_level=self.cpu_level
        )

        files = [f for f in glob.glob(os.path.join(self.input_dir, "*.fit*"))
                 if os.path.dirname(f) == self.input_dir]
        buckets = {}

        print("\n[~] Skimming headers and file sizes...")
        total_size_mb = 0
        valid_files = 0

        for f in files:
            if "Master_" in os.path.basename(f):
                continue
            meta = FITSHandler.extract_metadata(f, self.instrument)
            if not meta:
                continue

            key = f"{meta['object']}_{meta['filter']}_{meta['hardware_id']}_{meta['exposure']}s"
            # Epoch binning keeps separate nights separate. Without it, five
            # nights of the same field collapse into one master and the time
            # axis is averaged away before phase 3 ever sees it.
            epoch = epoch_key(meta, self.config.phase2.epoch_bin, self.logger)
            if epoch:
                key = f"{key}_{epoch}"
                meta["epoch"] = epoch
            buckets.setdefault(key, TargetGroup(key, meta)).raw_files.append(f)

            total_size_mb += os.path.getsize(f) / (1024 * 1024)
            valid_files += 1

        self.target_groups = list(buckets.values())
        if not self.target_groups:
            raise SystemExit("[ERROR] No valid FITS found.")

        max_stack = max(len(g.raw_files) for g in self.target_groups)
        self.hw.print_estimations(valid_files, max_stack, total_size_mb)

        self.run_dir = self.output_dir
        os.makedirs(self.run_dir, exist_ok=True)

        self.logger = get_logger("cassa_integrate", run_dir=self.run_dir)
        self.logger.info(f"Phase 2 output -> {self.run_dir}")
        self.logger.info("=" * 50)
        self.logger.info(f"PHASE 2: UNIFIED PIPELINE (Stacks + WCS + QA) - Cores: {cores}")
        self.logger.info(f"Instrument profile: {self.instrument.name}")
        self.logger.info("=" * 50)


    # --- Frame survey ---------------------------------------------------------
    def _survey_frames(self, paths, cfg):
        """What each frame is worth: its noise, sharpness, depth and conditions."""
        surveyed = []
        for path in paths:
            try:
                data = FITSHandler.load_data(path)
                subtracted, _ = MathEngine.extract_2d_background(
                    data, cfg.background_box, cfg.background_filter)
                noise = float(mad_std(subtracted, ignore_nan=True))
                header = fits.getheader(path)
            except Exception as exc:
                self.logger.warning("  [!] Could not read %s: %s",
                                    os.path.basename(path), exc)
                continue
            surveyed.append({
                "path": path,
                "noise": noise,
                "snr": float(np.nanpercentile(subtracted, 99) / noise) if noise > 0 else 0.0,
                # Measured in phase 1. Absent on frames reduced before that, and
                # every consumer here treats None as "unknown" rather than
                # substituting a guess.
                "fwhm": _positive_or_none(header.get("FWHMPX")),
                "exposure": _positive_or_none(header.get("EXPTIME")),
                "airmass": _positive_or_none(header.get("AIRMASS")),
            })
        return surveyed

    def _reject_poor_seeing(self, frames, cfg):
        """Drop frames far softer than the rest of the group.

        Nothing rejected frames on sharpness before: a badly defocused or
        wind-shaken frame was stacked in with full weight, broadening the master
        PSF for every frame that follows it.
        """
        measured = [fr["fwhm"] for fr in frames if fr["fwhm"]]
        if not cfg.fwhm_reject_factor or len(measured) < 3:
            return frames, []
        limit = float(np.median(measured)) * cfg.fwhm_reject_factor
        kept = [fr for fr in frames if not fr["fwhm"] or fr["fwhm"] <= limit]
        rejected = [fr for fr in frames if fr["fwhm"] and fr["fwhm"] > limit]
        return kept, rejected

    def _anchor_rank(self, frame, cfg):
        """Rank a frame as an anchor candidate: deep *and* sharp."""
        if not frame["fwhm"]:
            return frame["snr"]
        return frame["snr"] / (frame["fwhm"] ** cfg.anchor_sharpness_power)

    def _write_aligned_frames(self, group, stack, var_stack, dq_stack, exposures,
                              airmasses, backgrounds, fwhms, rejected, skipped,
                              wcs_solver):
        """Write each registered frame instead of one combined master.

        This is what ``phase2.steps.stack: false`` produces. The frames are
        aligned, flux-normalised and carry SCI/ERR/DQ exactly as they would
        going into the combine -- they simply are not averaged, so the time axis
        survives into phase 3. That is what a light curve needs: stacking a
        night of a variable star averages away the very signal being measured.

        Each product is named ``Aligned_<group>_NNNofM.fits`` and matches the
        glob phase 3 already uses (``*.fits``, excluding its own outputs), so
        the next phase needs no change to pick them up.

        Only the anchor is plate-solved. With registration on, every frame here
        shares the anchor's pixel grid *by construction*, so one solution is the
        solution for all of them and re-solving each would burn minutes per
        frame to rediscover the same WCS. With registration off that premise
        does not hold, and the frames are left unsolved rather than stamped with
        a WCS that is not theirs.
        """
        cfg = self.config.phase2
        total = len(stack)
        aligned_paths = []

        # strict=True: the three lists are appended to in lockstep, so a length
        # mismatch is a bug in this loop rather than something to truncate past.
        for index, (data, variance, dq) in enumerate(
                zip(stack, var_stack, dq_stack, strict=True), 1):
            error = np.sqrt(np.asarray(variance, dtype=float))
            error = np.where(np.isfinite(error), error, np.nan)

            stats = {
                "exposure_mean": exposures[index - 1],
                "airmass_mean": airmasses[index - 1],
                "background_level": backgrounds[0] if backgrounds else None,
                "fwhm": fwhms[index - 1],
                "n_fwhm_rejected": len(rejected),
                "epoch": group.meta.get("epoch"),
            }
            meta = dict(group.meta)
            meta["stack_stats"] = stats
            meta["steps_skipped"] = skipped

            path = os.path.join(
                self.run_dir,
                f"Aligned_{group.group_key}_{index:03d}of{total:03d}.fits",
            )
            # num_frames=1: each product is one frame, so TOT_EXP is its own
            # exposure and phase 3 divides by EXPMEAN as it always does.
            FITSHandler.save_master(
                np.asarray(data, dtype=float), error, _as_dq(dq, np.shape(data)),
                group.anchor_filepath, path, 1, meta,
            )
            aligned_paths.append(path)

        group.meta["aligned_paths"] = aligned_paths
        self.logger.info(
            "  [*] Stacking is switched off (phase2.steps); wrote %d aligned "
            "frame(s) instead of a master.", total,
        )

        if not cfg.steps.enabled("solve_wcs"):
            self.logger.info("  [*] WCS solving is switched off (phase2.steps).")
            return

        # Solve the anchor, then share its WCS with the frames registered onto it.
        group.master_filepath = aligned_paths[0]
        if not wcs_solver.solve(group, keep_temps=self.keep_temps):
            self.logger.warning("  [!] WCS solve failed for %s.",
                                os.path.basename(aligned_paths[0]))
            return
        if not cfg.steps.enabled("align"):
            self.logger.warning(
                "  [!] Alignment is off, so the anchor's WCS does not describe "
                "the other frames; they are left unsolved."
            )
            return
        self._share_wcs(aligned_paths[0], aligned_paths[1:])

    def _share_wcs(self, solved_path, targets):
        """Copy a solved WCS onto frames known to share its pixel grid."""
        from cassa_photometry.fits_utils import read_mef, write_mef

        try:
            _, _, _, solved_header = read_mef(solved_path)
        except Exception as exc:
            self.logger.warning("  [!] Could not read the solved WCS: %s", exc)
            return

        cards = [c for c in ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1",
                             "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                             "CDELT1", "CDELT2", "EQUINOX", "RADESYS", "ASTRMS")
                 if c in solved_header]
        if not cards:
            self.logger.warning("  [!] The solved frame carries no WCS cards to share.")
            return

        for path in targets:
            try:
                sci, err, dq, header = read_mef(path)
                for card in cards:
                    header[card] = solved_header[card]
                header["WCSFROM"] = (os.path.basename(solved_path),
                                     "WCS copied from this co-registered frame")
                write_mef(path, sci=sci, err=err, dq=dq, header=header)
            except Exception as exc:
                self.logger.warning("  [!] Could not share the WCS to %s: %s",
                                    os.path.basename(path), exc)
        self.logger.info("  [*] Shared the anchor's WCS with %d co-registered "
                         "frame(s).", len(targets))

    def _register_frame(self, t_data, t_err, t_dq, r_data, anchor, frame, cfg):
        """Bring one frame onto the anchor's pixel grid and flux scale.

        Returns ``(data, variance, dq, scale, n_stars, noise)``, or a leading
        ``None`` for a frame rejected on its flux scale.

        With ``phase2.steps.align`` off the frame is taken to be on that grid
        already -- the case for a well-tracked mount observing without dither.
        Two consequences are worth stating rather than discovering:

        * A shape mismatch is then an error for that frame, not something to
          resample away. Nothing else can be true if the grids are shared.
        * The transparency scale is measured from stars *matched between the
          two frames*, which is a by-product of solving for the transform. With
          no transform there are no matches, so the scale stays 1.0 and frames
          are combined un-normalised. In variable transparency that is a real
          loss of accuracy, which is why alignment is on by default.
        """
        if not cfg.steps.enabled("align"):
            if t_data.shape != r_data.shape:
                raise ValueError(
                    f"shape {t_data.shape} does not match the anchor's "
                    f"{r_data.shape}; frames must already share a pixel grid "
                    f"when phase2.steps.align is off"
                )
            noise = mad_std(t_data, ignore_nan=True)
            variance = _variance_from(t_err, t_data.shape, noise)
            return (t_data, variance, _as_dq(t_dq, t_data.shape), 1.0, 0, noise)

        tf, (src, _dst) = aa.find_transform(
            median_filter(t_data, 3), median_filter(r_data, 3),
            detection_sigma=cfg.align_detection_sigma, min_area=cfg.align_min_area,
        )
        # Size the aperture from the WORSE of the two frames: an aperture that
        # fits the sharp one loses flux from the soft one, and that loss is
        # exactly the seeing-dependent bias this measurement exists to avoid.
        worst_fwhm = max(
            anchor["fwhm"] or cfg.scale_aperture_fwhm_default,
            frame["fwhm"] or cfg.scale_aperture_fwhm_default,
        )
        scale = MathEngine.calc_scale(
            t_data, r_data, src, _dst,
            aperture_radius=worst_fwhm * cfg.scale_aperture_factor,
        )
        if not (cfg.scale_min <= scale <= cfg.scale_max):
            self.logger.warning(f"  [!] REJECTED: Scale {scale:.2f}x")
            return (None, None, None, scale, len(src), None)

        reg_data = MathEngine.register_plane(
            t_data, tf.inverse, r_data.shape, cfg.warp_order)
        norm_data = reg_data * scale
        noise = mad_std(norm_data, ignore_nan=True)

        t_var = _variance_from(t_err, t_data.shape, noise)
        reg_var = MathEngine.register_variance(t_var, tf.inverse, r_data.shape)
        norm_var = reg_var * (scale ** 2)

        # The DQ plane is resampled with the data, so a bad pixel stays a bad
        # pixel in the master's frame of reference.
        reg_dq = MathEngine.register_mask(
            _as_dq(t_dq, t_data.shape), tf.inverse, r_data.shape)

        return (norm_data, norm_var, reg_dq, scale, len(src), noise)

    def _remove_background(self, data, cfg):
        """Subtract the 2D sky, or keep it, per ``phase2.steps.subtract_background``.

        Returns ``(data, background)`` either way, so the caller records
        ``BKGLEVEL`` whichever branch ran -- the level is measured in both cases,
        it is only the *subtraction* that is optional.

        Switching this off is what extended-source work needs: the sky model
        follows large-scale structure, so subtracting it removes the outer disk
        of a resolved galaxy along with the sky. It is on by default because it
        is what makes registration and stacking well behaved.
        """
        subtracted, background = MathEngine.extract_2d_background(
            data, cfg.background_box, cfg.background_filter)
        if not cfg.steps.enabled("subtract_background"):
            return np.asarray(data, dtype=float), background
        return subtracted, background

    def _measure_master_fwhm(self, master):
        """The stack's own delivered PSF, measured on the stacked image.

        Returns None when ``phase2.steps.measure_fwhm`` is off, which is the
        same value a failed measurement yields -- so every consumer already
        handles it, and phase 3 falls back to ``phase3.fwhm``.
        """
        if not self.config.phase2.steps.enabled("measure_fwhm"):
            return None

        from cassa_photometry.phase4_diagnostics.psf import estimate_fwhm

        try:
            result = estimate_fwhm(
                np.asarray(master, dtype=float),
                fwhm_guess=self.config.phase4.fwhm_guess,
                threshold=self.config.phase4.detection_threshold,
                max_stars=self.config.phase4.max_stars,
                cutout=self.config.phase4.cutout,
            )
        except Exception:
            return None
        value = result.get("fwhm_px") if result else None
        return float(value) if value and np.isfinite(value) else None

    def _frame_weight(self, noise, fwhm, cfg):
        """Stacking weight for one frame.

        ``1/sigma**2`` maximises SNR for *extended* flux. A point source's SNR
        goes as ``1/(sigma*FWHM)`` because worse seeing spreads the same flux
        over more pixels, so the matching weight is ``1/(sigma*FWHM)**2``. Two
        frames of equal background noise but 1.5x different seeing currently
        contribute equally; they should not.

        Point-source weighting is the default because variable stars and
        supernovae are the science driver. Extended-source work wants
        ``stack_weight: extended``.
        """
        if not noise or noise <= 0:
            return 0.0
        weight = 1.0 / (noise ** 2)
        if cfg.stack_weight == "point_source" and fwhm:
            weight /= fwhm ** 2
        return weight

    def execute(self):
        cfg = self.config.phase2
        wcs_solver = WCSSolver(self.logger, self.config, instrument=self.instrument)
        failed_wcs = []

        # Say what will not be done, before doing anything. Phase 1 has said
        # this since the toggles were introduced; phase 2 stayed silent, which
        # is how five of its six toggles went years without a consumer.
        plan = resolve_plan("phase2", cfg.steps)
        names = [step.name for step in plan]
        # Whether the sky is removed before or after registration. True unless
        # the user explicitly ordered `align` first.
        bkg_first = not (
            "subtract_background" in names and "align" in names
            and names.index("align") < names.index("subtract_background")
        )
        self.logger.info("Phase 2 plan: %s", " -> ".join(names))

        skipped = cfg.steps.skipped()
        if skipped:
            self.logger.warning(
                "Deliberately skipping: %s. The masters record this as STEPSKIP "
                "so a partially-integrated stack cannot pass for a finished one.",
                ", ".join(skipped),
            )

        for group in self.target_groups:
            tot = len(group.raw_files)
            if tot < 2:
                continue

            self.logger.info(f"\n--- Stacking: {group.group_key} ({tot} frames) ---")

            # Anchor selection. The anchor sets both the astrometric reference
            # and the flux scale every other frame is normalised to, so it
            # should be the *best* frame -- which means sharp as well as deep.
            # Ranking on SNR alone lets a bright, soft frame win, and every
            # other frame is then matched to a blurred reference.
            surveyed = self._survey_frames(group.raw_files, cfg)
            if not surveyed:
                self.logger.warning("  [!] No readable frames in this group; skipping.")
                continue

            kept, rejected = self._reject_poor_seeing(surveyed, cfg)
            if rejected:
                self.logger.warning(
                    "  [!] Rejected %d frame(s) for poor seeing (FWHM > %.1fx the "
                    "group median): %s", len(rejected), cfg.fwhm_reject_factor,
                    ", ".join(os.path.basename(r["path"]) for r in rejected[:4]),
                )
            if len(kept) < 2:
                self.logger.warning("  [!] Fewer than 2 usable frames; skipping group.")
                continue

            anchor = max(kept, key=lambda fr: self._anchor_rank(fr, cfg))
            group.anchor_filepath = anchor["path"]
            best_noise = anchor["noise"]
            self.logger.info(
                "[+] Anchor set: %s (SNR %.0f, FWHM %s)",
                os.path.basename(group.anchor_filepath), anchor["snr"],
                f"{anchor['fwhm']:.2f} px" if anchor["fwhm"] else "unknown",
            )

            a_sci, a_err, a_dq = FITSHandler.load_planes(group.anchor_filepath)
            # The anchor is never warped, so for it the two orders are the same
            # operation; only the registration *reference* differs.
            r_data, r_bkg = self._remove_background(a_sci, cfg)
            r_var = _variance_from(a_err, r_data.shape, best_noise)

            stack = [r_data]
            var_stack = [r_var]
            bad_stack = [_bad_pixel_mask(a_dq, r_data.shape)]
            dq_stack = [_as_dq(a_dq, r_data.shape)]
            weights = [self._frame_weight(best_noise, anchor["fwhm"], cfg)]
            fwhms = [anchor["fwhm"]]
            exposures = [anchor["exposure"]]
            airmasses = [anchor["airmass"]]
            backgrounds = [float(np.nanmedian(r_bkg))]
            qa_stars, qa_scales, qa_noises = [], [], [best_noise]

            for frame in kept:
                f = frame["path"]
                if f == group.anchor_filepath:
                    continue
                t_sci, t_err, t_dq = FITSHandler.load_planes(f)
                # Resolved order decides whether the sky model is fitted to the
                # frame as observed or to the resampled one. Both are defensible
                # -- fitting before avoids interpolating the sky, fitting after
                # measures it on the grid the stack actually uses -- so the
                # registry links neither to the other and the user chooses.
                t_data = (self._remove_background(t_sci, cfg)[0] if bkg_first
                          else np.asarray(t_sci, dtype=float))

                try:
                    norm_data, norm_var, reg_dq, scale, n_stars, fn = \
                        self._register_frame(t_data, t_err, t_dq, r_data,
                                             anchor, frame, cfg)
                    if norm_data is None:
                        continue
                    if not bkg_first:
                        norm_data, _ = self._remove_background(norm_data, cfg)

                    stack.append(norm_data)
                    var_stack.append(norm_var)
                    dq_stack.append(reg_dq)
                    bad_stack.append(_bad_pixel_mask(reg_dq, r_data.shape))
                    weights.append(self._frame_weight(fn, frame["fwhm"], cfg))
                    fwhms.append(frame["fwhm"])
                    exposures.append(frame["exposure"])
                    airmasses.append(frame["airmass"])

                    qa_scales.append(scale)
                    qa_noises.append(fn)
                    if cfg.steps.enabled("align"):
                        qa_stars.append(n_stars)
                        self.logger.info(
                            "  [+] Aligned: %s | Stars: %d | Scale: %.2fx",
                            os.path.basename(f), n_stars, scale)
                    else:
                        self.logger.info("  [+] Added unregistered: %s",
                                         os.path.basename(f))
                except Exception as exc:
                    self.logger.error(f"  [-] Alignment failed for {os.path.basename(f)}: {exc}")

            group.successful_frames = len(stack)

            # Stacking off: write each registered frame instead of combining
            # them. `min_frames_to_stack` is a *combine* threshold -- it exists
            # because a handful of frames carries too little information to
            # reject an outlier -- so it must not gate a path that rejects
            # nothing. A single frame is a legitimate result here.
            if not cfg.steps.enabled("stack"):
                self._write_aligned_frames(
                    group, stack, var_stack, dq_stack, exposures, airmasses,
                    backgrounds, fwhms, rejected, skipped, wcs_solver,
                )
                continue

            if group.successful_frames < cfg.min_frames_to_stack:
                self.logger.warning(
                    "  [!] Only %d frame(s) survived alignment (need %d); "
                    "skipping this group.",
                    group.successful_frames, cfg.min_frames_to_stack,
                )
                continue

            self.logger.info("  [*] Integrating cube (inverse-variance weighted)...")
            cube = np.array(stack)
            varcube = np.array(var_stack)
            w_arr = np.array(weights) / sum(weights)
            use_clip = group.successful_frames >= cfg.sigma_clip_min_frames

            master, master_var, n_used = MathEngine.weighted_stack(
                cube, varcube, w_arr,
                sigma=cfg.stack_sigma, maxiters=cfg.stack_maxiters, use_sigma_clip=use_clip,
                bad_pixels=np.array(bad_stack) if cfg.mask_bad_pixels else None,
            )
            master_err = np.sqrt(master_var)
            # An uncovered pixel has unknown, not zero, uncertainty. Writing 0
            # gives it infinite inverse-variance weight everywhere downstream.
            master_err = np.where(np.isfinite(master_err), master_err, np.nan)

            # Carry phase 1's data quality into the master. Rebuilding DQ from
            # scratch here -- which is what this used to do -- discards every
            # saturation, cosmic-ray and bad-pixel flag the reduction found, so
            # nothing downstream can exclude a saturated star from a zero point.
            master_dq = _combine_dq(dq_stack, n_used, master)
            group.meta["dq_counts"] = _dq_summary(master_dq)

            # What the stack was made of, so the master states its own
            # provenance instead of leaving phase 3 to guess.
            measured_fwhm = self._measure_master_fwhm(master)
            group.meta["stack_stats"] = {
                "exposure_mean": _mean_of(exposures),
                "airmass_mean": _mean_of(airmasses),
                "background_level": _mean_of(backgrounds),
                # Measured ON the master, not averaged from the frames:
                # alignment residuals broaden a stack, so the combined PSF is
                # worse than the mean of its inputs and only direct measurement
                # gives the number phase 3 actually needs.
                "fwhm": measured_fwhm,
                "fwhm_std": _std_of(fwhms),
                "fwhm_min": _min_of(fwhms),
                "fwhm_max": _max_of(fwhms),
                "n_fwhm_rejected": len(rejected),
                "epoch": group.meta.get("epoch"),
            }
            group.meta["steps_skipped"] = skipped
            group.meta["steps_plan"] = names

            group.master_filepath = os.path.join(
                self.run_dir, f"Master_{group.group_key}_{group.successful_frames}fr.fits"
            )
            FITSHandler.save_master(
                master, master_err, master_dq,
                group.anchor_filepath, group.master_filepath,
                group.successful_frames, group.meta,
            )

            self._log_telemetry(group, r_data, master, qa_stars, qa_scales, qa_noises)
            del stack, var_stack, cube, varcube, master, master_var, master_err
            gc.collect()

            if not cfg.steps.enabled("solve_wcs"):
                self.logger.info("  [*] WCS solving is switched off (phase2.steps).")
                continue

            self.logger.info("  [*] Initiating Astrometry WCS Solve...")
            if not wcs_solver.solve(group, keep_temps=self.keep_temps):
                failed_wcs.append(group)

        if failed_wcs and wcs_solver.global_anchor_ra is not None:
            self.logger.info("\n--- INITIATING WCS RESCUE PASS ---")
            for g in failed_wcs:
                wcs_solver.solve(g, is_rescue=True, keep_temps=self.keep_temps)

        if cfg.steps.enabled("visual_qa"):
            VisualQAGenerator.generate_pdf(self.target_groups, self.run_dir,
                                           self.logger, self.config)
        else:
            self.logger.info("  [*] Visual QA PDF is switched off (phase2.steps).")

        self.logger.info("\n" + "=" * 50)
        self.logger.info("PHASE 2 PIPELINE COMPLETE.")

    def _log_telemetry(self, group, r_data, master, qa_stars, qa_scales, qa_noises):
        """Emit the automated QA telemetry (SNR boost / integration efficiency)."""
        self.logger.info("\n   >>> AUTOMATED QA TELEMETRY <<<")
        if qa_stars:
            self.logger.info(f"   [*] STAR MATCH : Avg: {int(np.mean(qa_stars))} | Min: {int(np.min(qa_stars))}")
        if qa_scales:
            self.logger.info(f"   [*] FLUX SCALE : Avg: {np.mean(qa_scales):.3f}x")

        c_anchor = MathEngine.center_crop(r_data)
        c_master = MathEngine.center_crop(master)
        anchor_noise = mad_std(c_anchor, ignore_nan=True)
        master_noise = mad_std(c_master, ignore_nan=True)

        actual_boost = anchor_noise / master_noise if master_noise > 0 else 0
        variances = np.array(qa_noises) ** 2
        expected_noise = np.sqrt(1.0 / np.sum(1.0 / variances))
        weighted_max = anchor_noise / expected_noise if expected_noise > 0 else 0
        efficiency = (actual_boost / weighted_max) * 100 if weighted_max > 0 else 0

        self.logger.info(f"   [*] PERFORMANCE: SNR Boost: {actual_boost:.2f}x | Efficiency: {efficiency:.1f}%")
        if efficiency < 50.0:
            self.logger.warning("       -> WARNING: Integration efficiency critically low (<50%).")
        self.logger.info("   --------------------------------------\n")


def run(input_dir, output_dir=None, keep_temps=False, config=None, logger=None,
        assume_yes=False, instrument=None):
    """Set up and run phase 2 over a directory of calibrated frames."""
    input_dir = os.path.abspath(input_dir)
    if not os.path.exists(input_dir):
        raise SystemExit(f"Directory {input_dir} does not exist.")

    pipeline = IntegrationPipeline(input_dir, output_dir=output_dir, keep_temps=keep_temps,
                                   config=config, logger=logger, instrument=instrument,
                                   assume_yes=assume_yes)
    pipeline.setup()

    if not assume_yes and pipeline.hw._interactive():
        try:
            proceed = input("\nProceed with Unified Integration? [Y/n]: ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            proceed = ""
        if proceed not in ("y", ""):
            return
    pipeline.execute()
