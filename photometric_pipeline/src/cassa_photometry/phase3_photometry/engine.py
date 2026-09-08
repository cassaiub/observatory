"""Photometry engine: zero point, flux calibration, and error-carrying catalogs.

Every measurement now carries an uncertainty. The zero point is a sigma-clipped,
inverse-variance weighted combination of per-star offsets and comes with a
``MAGZERR``; the output catalog exposes ``FLUXERR_ISO``, ``MAGERR_ISO``, ``SNR``
and DQ ``FLAGS``. The ``ERR`` plane written by phases 1-2 drives the photometric
errors; if it is missing (legacy input) it is reconstructed from the background.
"""

import astropy.units as u
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import SigmaClip, sigma_clip, sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from photutils.segmentation import SourceCatalog, deblend_sources, detect_sources
from photutils.utils import calc_total_error
from scipy import ndimage

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import read_mef, write_mef
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.phase3_photometry import apcor, morphology, photsys
from cassa_photometry.phase3_photometry.catalogs import fetch_reference_catalog

# 2.5 / ln(10): converts a fractional flux error into a magnitude error.
_POGSON = 2.5 / np.log(10.0)


class UniversalPhotometryEngine:
    def __init__(self, fwhm_estimate=None, detection_threshold=None, config=None, logger=None):
        self.config = config or load_config()
        cfg = self.config.phase3
        # An explicit --fwhm overrides the measurement; otherwise the frame's
        # own FWHMPX is used and this is only the last-resort fallback.
        self._fwhm_override = fwhm_estimate
        self.fwhm = fwhm_estimate if fwhm_estimate is not None else cfg.fwhm
        self.threshold = detection_threshold if detection_threshold is not None else cfg.detection_threshold
        self.logger = logger or get_logger("cassa_photometry")
        self.zero_point = None
        self.zero_point_err = None
        self.n_zp_stars = 0
        #: The zero point before the aperture correction, kept for the header.
        self.aperture_zero_point = None
        self.aperture_correction = apcor.ApertureCorrection(0.0, method="unmeasured")
        #: Single-frame exposure the instrumental magnitudes are normalised by.
        self.exposure = 1.0
        self.band = None
        self._classification = None

    # -- helpers --------------------------------------------------------------
    def _load(self, path):
        """Return ``(sci, err, dq, header)``; reconstruct ERR if it is missing."""
        sci, err, dq, header = read_mef(path)
        if err is None:
            _, _, std = sigma_clipped_stats(sci, sigma=3.0)
            # A master is a weighted MEAN of N frames, so its effective gain is
            # about N: each electron of the mean was measured N times. Assuming
            # 1 overestimates the Poisson term by roughly sqrt(N), and this
            # fallback feeds every magnitude error and every zero-point weight.
            n_frames = header.get("STACKCNT")
            try:
                effective_gain = max(float(n_frames), 1.0)
            except (TypeError, ValueError):
                effective_gain = 1.0
            err = calc_total_error(np.clip(sci, 0, None), np.full_like(sci, std),
                                   effective_gain=effective_gain)
            self.logger.warning(
                "No ERR plane found; reconstructed background+Poisson errors "
                "with an effective gain of %.3g.", effective_gain,
            )
        return sci, err, dq, header

    # -- zero point -----------------------------------------------------------
    def resolve_fwhm(self, header):
        """The FWHM to size apertures from, in pixels.

        Read from the frame rather than configured. ``FWHMPX`` is measured by
        phases 1 and 2 from the stars in the image; falling back to a config
        default means a **fixed aperture whatever the seeing was**, which is a
        per-epoch error in the zero point of up to 0.3 mag -- precisely the
        thing a light curve cannot tolerate.

        Note that ``SEEING``, if present, is *not* used: it is the DIMM's
        zenith-corrected atmospheric seeing at a reference wavelength, and the
        delivered PSF is larger by the airmass, wavelength, guiding and optics
        terms. Sizing apertures from it under-sizes them by roughly 30%.
        """
        if self._fwhm_override is not None:
            return float(self._fwhm_override)
        measured = header.get("FWHMPX")
        try:
            value = float(measured)
            if np.isfinite(value) and value > 0:
                return value
        except (TypeError, ValueError):
            pass
        self.logger.warning(
            "No FWHMPX card; falling back to phase3.fwhm = %.2f px. Apertures "
            "will not be matched to this frame's seeing.", self.config.phase3.fwhm,
        )
        return float(self.config.phase3.fwhm)

    def search_radius(self, wcs, shape):
        """Radius of the reference-catalog cone, covering the whole field.

        Two corrections over ``nx/2 * CD1_1``: ``proj_plane_pixel_scales``
        instead of ``CD1_1`` (which is ``scale*cos(theta)`` on a rotated field),
        and the **half-diagonal** instead of the half-width. Measured on a real
        master, the old expression asked for 5.04' where 7.13' was needed --
        half the field area had no reference stars available at all.
        """
        ny, nx = shape
        scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial)))  # deg/pixel
        return (0.5 * np.hypot(nx, ny) * scale) * u.deg

    def calculate_local_zero_point(self, master_science_fits, science_band="R"):
        """Compute a filter-wise zero point with uncertainty (differential photometry)."""
        cfg = self.config.phase3
        self.logger.info("Calculating local zero point (differential photometry)...")
        data, err, dq, header = self._load(master_science_fits)
        wcs = WCS(header)
        ny, nx = data.shape

        self.fwhm = self.resolve_fwhm(header)
        self.exposure = _single_frame_exposure(header, self.logger)
        self.band = str(science_band).upper()

        center = wcs.pixel_to_world(nx / 2, ny / 2)
        catalog = fetch_reference_catalog(
            center, self.search_radius(wcs, data.shape), science_band, self.logger,
            config=self.config,
        )
        cat_coords = SkyCoord(ra=catalog["RAJ2000"], dec=catalog["DEJ2000"],
                              unit=(u.deg, u.deg))

        # One background model, shared with the catalog path. Three different
        # notions of "sky" in one module is how two of them end up wrong.
        background = _background_of(data, self.config)
        subtracted = data - background.background
        threshold = self.threshold * background.background_rms

        finder = DAOStarFinder(fwhm=self.fwhm, threshold=float(np.nanmedian(threshold)))
        sources = finder.find_stars(subtracted)
        if sources is None or len(sources) == 0:
            raise ValueError("No stars detected to determine a zero point.")

        positions = np.column_stack([sources["xcentroid"], sources["ycentroid"]])
        det_coords = wcs.pixel_to_world(positions[:, 0], positions[:, 1])

        pairs = _mutual_matches(det_coords, cat_coords, self._match_radius(header))
        if not pairs:
            raise ValueError("Could not match any stars to the catalog for a zero point.")
        det_index = np.array([p[0] for p in pairs])
        cat_index = np.array([p[1] for p in pairs])

        flux, flux_err, usable = self._aperture_flux(
            subtracted, err, dq, positions[det_index]
        )
        if not np.any(usable):
            raise ValueError("No usable calibrator stars survived the quality cuts.")

        cat_mag = np.asarray(catalog["Ref_Mag"], dtype=float)[cat_index][usable]
        cat_mag_err = np.asarray(catalog["Ref_Mag_Err"], dtype=float)[cat_index][usable]

        # Instrumental magnitude per SECOND, so the zero point does not silently
        # absorb the exposure time and stays comparable between epochs.
        inst_mag = -2.5 * np.log10(flux[usable] / self.exposure)
        inst_mag_err = _POGSON * flux_err[usable] / flux[usable]

        zps = cat_mag - inst_mag
        zp_errs = np.hypot(inst_mag_err, cat_mag_err)

        aperture_zp, aperture_zp_err, self.n_zp_stars = combine_zeropoints(
            zps, zp_errs, sigma=cfg.zp_sigma_clip,
        )

        # From the aperture's flux scale to a total-flux one. Without this the
        # zero point and the catalog magnitudes it calibrates are measured on
        # different scales, and the error is a tilt rather than an offset.
        self.aperture_correction = self._measure_aperture_correction(
            subtracted, positions[det_index][usable], data.shape
        )
        self.zero_point = aperture_zp + self.aperture_correction.value
        # A scalar correction over a field-dependent PSF is right at the centre
        # and wrong in the corners; the spread is an honest part of the error.
        extra = (self.aperture_correction.scatter
                 if self.aperture_correction.method != "surface" else 0.0)
        self.zero_point_err = float(np.hypot(aperture_zp_err, extra))
        self.aperture_zero_point = aperture_zp

        self.logger.info(
            "Zero point (%s-band): %.4f +/- %.4f mag from %d stars "
            "[aperture %.4f, apcor %+.4f, %s].",
            science_band, self.zero_point, self.zero_point_err, self.n_zp_stars,
            aperture_zp, self.aperture_correction.value, photsys.describe(self.band),
        )
        return self.zero_point

    def _match_radius(self, header):
        """Cross-match radius, derived from the astrometric fit where possible."""
        from cassa_photometry.phase3_photometry.verify import match_radius_from_header

        radius, why = match_radius_from_header(header, self.config)
        self.logger.info("    Cross-match radius: %.2f\" [%s]", radius, why)
        return float(radius) * u.arcsec

    def _aperture_flux(self, data, err, dq, positions):
        """Background-subtracted aperture flux for every calibrator, vectorised.

        Three changes from the original, all of which bias a zero point:

        * the annulus background is a **sigma-clipped median**, not a mean, so a
          neighbouring star or an un-flagged cosmic ray in the annulus does not
          bias every calibrator it touches;
        * the annulus is moved **outside the PSF wings** -- at 3-4x FWHM a
          Moffat still puts about 1.25% of the star's own light there, which was
          being subtracted off as sky;
        * stars with a **bad DQ pixel in the aperture are rejected**. The
          brightest catalog stars saturate first and are exactly the ones
          inverse-variance weighting trusts most.
        """
        from photutils.aperture import ApertureStats

        cfg = self.config.phase3
        aperture = CircularAperture(positions, r=self.fwhm * cfg.aperture_r_factor)
        annulus = CircularAnnulus(positions,
                                  r_in=self.fwhm * cfg.annulus_in_factor,
                                  r_out=self.fwhm * cfg.annulus_out_factor)

        photometry = aperture_photometry(data, aperture, error=err)
        raw = np.asarray(photometry["aperture_sum"], dtype=float)
        raw_err = (np.asarray(photometry["aperture_sum_err"], dtype=float)
                   if "aperture_sum_err" in photometry.colnames
                   else np.zeros(len(raw)))

        stats = ApertureStats(data, annulus, sigma_clip=SigmaClip(sigma=3.0, maxiters=5))
        sky = np.asarray(stats.median, dtype=float)
        sky_std = np.asarray(stats.std, dtype=float)
        n_sky = np.asarray(stats.sum_aper_area.value if hasattr(stats.sum_aper_area, "value")
                           else stats.sum_aper_area, dtype=float)

        flux = raw - sky * aperture.area
        # The subtracted sky is itself uncertain; a median of n pixels has a
        # standard error of about 1.25*sigma/sqrt(n).
        with np.errstate(invalid="ignore", divide="ignore"):
            sky_err = 1.2533 * sky_std / np.sqrt(np.clip(n_sky, 1, None))
        flux_err = np.sqrt(raw_err**2 + (aperture.area * sky_err) ** 2)

        usable = np.isfinite(flux) & (flux > 0) & np.isfinite(flux_err) & (flux_err > 0)
        if dq is not None and cfg.zp_reject_flagged:
            usable &= ~_aperture_has_bad_pixels(dq, aperture)
        rejected = int(np.count_nonzero(~usable))
        if rejected:
            self.logger.info(
                "    %d of %d calibrator(s) rejected (saturated, flagged or "
                "non-positive flux).", rejected, len(usable),
            )
        return flux, flux_err, usable

    def _measure_aperture_correction(self, data, positions, shape):
        """Curve of growth on the calibrators, as a surface where possible."""
        cfg = self.config.phase3
        if not cfg.steps.enabled("aperture_correction"):
            return apcor.ApertureCorrection(0.0, method="disabled")
        return apcor.measure(
            data, positions, self.fwhm,
            aperture_radius=self.fwhm * cfg.aperture_r_factor,
            shape=shape, plateau_factor=cfg.apcor_total_factor, logger=self.logger,
        )

    # -- flux-calibrated image ------------------------------------------------
    def export_flux_calibrated_image(self, target_fits, output_filename):
        """Write a SCI/ERR image in Jy/pixel (propagating the error plane).

        Two corrections over the original:

        * **The AB offset is applied.** The conversion is the AB relation, but
          B/V/R/I are calibrated against Vega-based magnitudes. Feeding a Vega
          magnitude through it unchanged is wrong by the band's AB offset --
          about 8% in B.
        * **The unit is ``Jy/pixel``, not ``Jy``.** The scaling is per pixel; a
          source's flux is the sum over its pixels. Labelling a per-pixel value
          ``Jy`` puts anyone doing surface photometry off by the pixel count.
        """
        self._require_zp()
        cfg = self.config.phase3
        data, err, dq, header = self._load(target_fits)

        # The zero point calibrates magnitudes in the band's own system; the Jy
        # conversion needs AB.
        offset = photsys.ab_offset(self.band) if cfg.apply_ab_offset else 0.0
        ab_zero_point = self.zero_point + offset
        conversion = cfg.ab_flux_zero_jy * 10 ** (-ab_zero_point / 2.5)
        self.logger.info(f"Writing flux-calibrated image: {output_filename}")

        header["BUNIT"] = "Jy/pixel"
        header["MAGZERO"] = (self.zero_point, "Total-flux zero point [mag]")
        header["MAGZERR"] = (self.zero_point_err, "Zero point uncertainty [mag]")
        header["NZPSTARS"] = (self.n_zp_stars, "Number of stars used for zero point")
        header["FLUXCAL"] = (conversion, "Jy/pixel per data unit")
        self._stamp_calibration(header, offset, ab_zero_point)

        write_mef(output_filename, sci=data * conversion,
                  err=(None if err is None else err * conversion), dq=dq, header=header,
                  history=f"Flux calibrated to Jy/pixel using ZP {ab_zero_point:.4f} (AB)")

    def _stamp_calibration(self, header, ab_offset_mag, ab_zero_point):
        """Record what system the numbers are on, and how they got there.

        A product that does not state its own calibration cannot be combined
        with another safely, and mixed vintages are then undetectable.
        """
        from cassa_photometry.fits_utils import CALVERS

        header["CALVERS"] = (CALVERS, "Calibration vintage")
        header["PHOTSYS"] = (photsys.system_of(self.band),
                             "System the calibrated magnitudes are on")
        # What ran, and what did not. The catalog is a CSV and cannot carry
        # this, so the flux-calibrated image is where phase 3's plan is on
        # record -- the same contract phase 1 keeps with CALSKIP.
        from cassa_photometry.steps import resolved_names

        steps = self.config.phase3.steps
        header["STEPPLAN"] = (",".join(resolved_names("phase3", steps)),
                              "Phase 3 steps run, in order")
        if steps.skipped():
            header["STEPSKIP"] = (",".join(steps.skipped()),
                                  "Phase 3 steps deliberately skipped")
        header["ABOFFSET"] = (round(float(ab_offset_mag), 4),
                              "[mag] AB minus Vega, applied to the Jy conversion")
        header["ZPAB"] = (round(float(ab_zero_point), 4), "Zero point on the AB system")
        header["APCOR"] = (round(float(self.aperture_correction.value), 4),
                           "[mag] Aperture -> total flux correction")
        header["APCORMTH"] = (self.aperture_correction.method,
                              "How APCOR was determined")
        header["APCORRMS"] = (round(float(self.aperture_correction.scatter), 4),
                              "[mag] APCOR spread across the field")
        if self.aperture_zero_point is not None:
            header["ZPAPER"] = (round(float(self.aperture_zero_point), 4),
                                "Zero point before the aperture correction")
        header["EXPNORM"] = (round(float(self.exposure), 4),
                             "[s] Exposure the instrumental mags were divided by")
        header["FWHMUSED"] = (round(float(self.fwhm), 4),
                              "[pixel] FWHM the apertures were sized from")
        wavelength = photsys.effective_wavelength_nm(self.band)
        if wavelength:
            header["PHOTWAVE"] = (wavelength, "[nm] Effective wavelength of the band")
        # A Jy value is a monochromatic flux density at that wavelength, valid
        # for a source with the calibrators' spectral shape. Stating it beats
        # implying it.
        header["PHOTREF"] = ("calibrator SED", "Jy conversion assumes this spectral shape")
        # The master sits on the anchor frame's throughput scale, not an
        # absolute one, which is exactly why the zero point must stay per-epoch.
        header["ZPSCOPE"] = ("per-epoch", "ZP is tied to this stack; do not average")

    # -- source catalog -------------------------------------------------------
    def generate_full_catalog(self, target_fits, output_csv, output_segmap=None):
        """Detect sources and write a catalog with flux/magnitude uncertainties.

        A zero point is *not* required. Without one -- a narrowband filter has no
        broadband reference catalog to calibrate against -- detection, fluxes and
        instrumental magnitudes are still valid and still written;
        ``MAG_ISO`` comes back NaN rather than carrying a fabricated
        calibration.
        """
        cfg = self.config.phase3
        self.logger.info("Generating full object catalog and segmentation map...")
        data, err, dq, header = self._load(target_fits)
        wcs = WCS(header)

        self.fwhm = self.resolve_fwhm(header)
        self.exposure = _single_frame_exposure(header, self.logger)

        bkg = _background_of(data, self.config)
        data_sub = data - bkg.background
        threshold = self.threshold * bkg.background_rms

        kernel = Gaussian2DKernel(x_stddev=self.fwhm / 2.35)
        kernel.normalize()
        convolved = convolve(data_sub, kernel)

        segmap = detect_sources(convolved, threshold, npixels=cfg.detect_npixels)
        if segmap is None:
            self.logger.warning("No sources detected.")
            return
        # Deblend on the same (convolved) data the detection used. Deblending
        # the unconvolved image splits sources on noise peaks that the detection
        # never saw.
        segmap = deblend_sources(convolved, segmap, npixels=cfg.detect_npixels,
                                 nlevels=cfg.deblend_nlevels, contrast=cfg.deblend_contrast)
        self.logger.info(f"Detected {segmap.nlabels} objects after deblending.")

        if output_segmap:
            seg_header = header.copy()
            seg_header["BUNIT"] = "ID"
            write_mef(output_segmap, sci=segmap.data.astype(np.int32), header=seg_header,
                      history="Segmentation map from UniversalPhotometryEngine")

        cat = SourceCatalog(data_sub, segmap, error=err, wcs=wcs,
                            convolved_data=convolved)
        tbl = cat.to_table()

        flux = np.asarray(tbl["segment_flux"], dtype=float)
        flux_err = np.asarray(tbl["segment_fluxerr"], dtype=float)

        # Keep non-positive fluxes. Dropping them -- which is what this used to
        # do -- biases faint number counts, distorts the limiting magnitude, and
        # makes upper limits impossible: a non-detection is a measurement, and
        # for time-domain work it is the measurement that constrains a rise.
        measurable = np.isfinite(flux) & np.isfinite(flux_err) & (flux_err > 0)
        if not cfg.keep_negative_flux:
            measurable &= flux > 0
        n_dropped = int(np.count_nonzero(~measurable))
        if n_dropped:
            self.logger.info("    %d source(s) had no usable measurement.", n_dropped)
        tbl, cat = tbl[measurable], cat[measurable]
        flux, flux_err = flux[measurable], flux_err[measurable]

        positive = flux > 0
        with np.errstate(invalid="ignore", divide="ignore"):
            inst_mag = np.where(positive, -2.5 * np.log10(np.abs(flux) / self.exposure), np.nan)
            inst_mag_err = np.where(positive, _POGSON * flux_err / np.abs(flux), np.nan)
        # A 3-sigma upper limit, defined for every source including the
        # non-detections, which is the point of keeping them.
        limit_mag = -2.5 * np.log10(3.0 * flux_err / self.exposure)

        if self.zero_point is None:
            mag = np.full_like(inst_mag, np.nan)
            mag_err = np.full_like(inst_mag, np.nan)
            limit_mag = np.full_like(limit_mag, np.nan)
        else:
            # With a fitted surface the correction is applied per source, since
            # the PSF -- and so the fraction an aperture catches -- varies across
            # the field. With a scalar, it is already inside self.zero_point.
            if self.aperture_correction.method == "surface":
                base = self.aperture_zero_point if self.aperture_zero_point is not None \
                    else self.zero_point
                correction = self.aperture_correction.at(
                    np.asarray(tbl["xcentroid"], dtype=float),
                    np.asarray(tbl["ycentroid"], dtype=float),
                )
                mag = inst_mag + base + np.asarray(correction, dtype=float)
            else:
                mag = inst_mag + self.zero_point
            # MAGERR_STAT is the random part; MAGZERR is a systematic shared by
            # every source in the frame, so adding it per-source (as MAGERR_ISO
            # does, for compatibility) double-counts it in any later average.
            mag_err = np.hypot(inst_mag_err, self.zero_point_err)
            limit_mag = limit_mag + self.zero_point

        # Column names follow the SExtractor/SEP convention, which is what
        # anyone receiving this catalog will expect. Two consequences are
        # deliberate: the measured quantity is isophotal, so it is named ISO
        # rather than left to be mistaken for a total magnitude; and X/Y_IMAGE
        # are 1-indexed per the FITS convention, unlike photutils' 0-indexed
        # centroids, so the positions overlay correctly in DS9 and friends.
        # A magnitude measured in the SAME aperture the zero point was, so at
        # least one catalog magnitude is exactly consistent with the calibration.
        # MAG_ISO is isophotal -- its name says so -- and is not a total
        # magnitude; the two differ by a brightness-dependent amount.
        (aperture_flux, aperture_flux_err, aperture_mag,
         aperture_mag_err) = self._aperture_magnitudes(data_sub, err, tbl)

        out = Table()
        out["NUMBER"] = tbl["label"]
        out["ALPHA_J2000"] = cat.sky_centroid.ra.deg
        out["DELTA_J2000"] = cat.sky_centroid.dec.deg
        out["X_IMAGE"] = np.asarray(tbl["xcentroid"], dtype=float) + 1.0
        out["Y_IMAGE"] = np.asarray(tbl["ycentroid"], dtype=float) + 1.0
        out["FLUX_ISO"] = flux
        out["FLUXERR_ISO"] = flux_err
        out["MAG_INST"] = inst_mag
        out["MAGERR_INST"] = inst_mag_err
        out["MAG_ISO"] = mag
        out["MAGERR_ISO"] = mag_err
        out["MAGERR_STAT"] = inst_mag_err
        out["MAG_APER"] = aperture_mag
        out["MAGERR_APER"] = aperture_mag_err
        out["FLUX_APER"] = aperture_flux
        out["FLUXERR_APER"] = aperture_flux_err
        out["LIMIT_MAG"] = limit_mag
        out["SNR"] = np.where(flux_err > 0, flux / flux_err, np.nan)
        out["ISOAREA_IMAGE"] = tbl["area"]

        out["FLAGS"] = _segment_flags(segmap.data, np.asarray(tbl["label"]), dq)

        # Kron ("AUTO") is the standard total magnitude; PSF is the matched
        # filter for a point source. Their difference is what classifies.
        kron_flux, kron_err, kron_mag, kron_mag_err, kron_radius = self._kron_magnitudes(cat)
        psf_mag, psf_mag_err, psf_chi2 = self._psf_magnitudes(data_sub, err, tbl, header)

        out["MAG_AUTO"] = kron_mag
        out["MAGERR_AUTO"] = kron_mag_err
        out["FLUX_AUTO"] = kron_flux
        out["FLUXERR_AUTO"] = kron_err
        out["KRON_RADIUS"] = kron_radius
        out["MAG_PSF"] = psf_mag
        out["MAGERR_PSF"] = psf_mag_err
        out["CHI2_PSF"] = psf_chi2
        out["PSF_MINUS_AUTO"] = psf_mag - kron_mag
        out["FLUX_RADIUS"] = _flux_radius(cat)
        out["FWHM_IMAGE"] = _source_fwhm(tbl, cat)

        if self.config.phase3.steps.enabled("classification"):
            classification = morphology.classify(
                out["PSF_MINUS_AUTO"], out["MAG_ISO"] if self.zero_point is not None
                else out["MAG_INST"],
                np.asarray(out["SNR"], dtype=float),
                saturated=_has_flag(out["FLAGS"], "saturated"),
                edge=_near_edge(tbl, data.shape, self.fwhm),
                logger=self.logger,
            )
            out["CLASS"] = classification.classes
            out["CLASS_STAR"] = classification.stellarity
            # The best available total magnitude: PSF for point sources, Kron
            # for extended ones. This is the column to use once it exists.
            out["MAG_BEST"] = np.where(classification.classes == morphology.STAR,
                                       psf_mag, kron_mag)
            out["MAGERR_BEST"] = np.where(classification.classes == morphology.STAR,
                                          psf_mag_err, kron_mag_err)
        else:
            # Classification off. CLASS and CLASS_STAR are omitted rather than
            # filled with a placeholder: a column of "UNKNOWN" invites being
            # read as a measurement, and any downstream cut on it would then be
            # silently meaningless.
            #
            # MAG_BEST cannot be omitted -- it is the documented column to use
            # -- so it falls back to Kron for every source. Kron is the general
            # total magnitude and is correct for extended sources; on a point
            # source it is noisier than the PSF fit rather than biased, which
            # is the safe direction to err.
            #
            # The catalog is a CSV and carries no header cards, so the *absence*
            # of CLASS/CLASS_STAR is what records this: a consumer that cuts on
            # either finds the column missing and fails loudly, rather than
            # cutting on a placeholder and silently selecting nothing.
            classification = None
            out["MAG_BEST"] = kron_mag
            out["MAGERR_BEST"] = kron_mag_err
            self.logger.info(
                "Star/galaxy classification is switched off (phase3.steps); "
                "CLASS/CLASS_STAR omitted and MAG_BEST falls back to Kron "
                "(MAG_AUTO) for every source."
            )
        # ELLIPTICITY stays, demoted from decision-maker to one feature among
        # several: it is shape, not concentration.
        out["ELLIPTICITY"] = _ellipticity(tbl, cat)
        self._classification = classification

        if self.zero_point is not None:
            # Which column to actually use. MAG_ISO is isophotal: it captures a
            # brightness-dependent fraction of a source, so an aperture-derived
            # zero point applied to it produces a TILT, not an offset. Measured
            # against simulated truth, MAG_ISO runs from +0.18 mag at V=12 to
            # +1.70 at V=16.5 (a slope of 0.35 mag per mag), while MAG_APER is
            # flat to -0.007 +/- 0.05 mag over the same range.
            self.logger.info(
                "    Use MAG_APER for photometry. MAG_ISO is isophotal and is "
                "not a total magnitude; it is retained for compatibility."
            )

        out.write(output_csv, format="csv", overwrite=True)
        calibration = ("instrumental magnitudes only -- no zero point"
                       if self.zero_point is None else "calibrated")
        self.logger.info(f"Catalog saved to {output_csv} ({len(out)} sources, {calibration}).")

    def _aperture_magnitudes(self, data, err, tbl):
        """Flux and magnitude in the zero point's own aperture.

        This is the column that makes the catalog and the calibration agree.
        The isophotal ``segment_flux`` the zero point used to be applied to
        captures a *brightness-dependent* fraction of a source, so applying an
        aperture zero point to it tilts the magnitude scale.
        """
        from photutils.aperture import ApertureStats

        cfg = self.config.phase3
        positions = np.column_stack([
            np.asarray(tbl["xcentroid"], dtype=float),
            np.asarray(tbl["ycentroid"], dtype=float),
        ])
        empty = np.full(len(positions), np.nan)
        if len(positions) == 0:
            return empty, empty, empty, empty

        aperture = CircularAperture(positions, r=self.fwhm * cfg.aperture_r_factor)
        annulus = CircularAnnulus(positions,
                                  r_in=self.fwhm * cfg.annulus_in_factor,
                                  r_out=self.fwhm * cfg.annulus_out_factor)
        try:
            photometry = aperture_photometry(np.nan_to_num(data, nan=0.0), aperture,
                                             error=err)
            raw = np.asarray(photometry["aperture_sum"], dtype=float)
            raw_err = (np.asarray(photometry["aperture_sum_err"], dtype=float)
                       if "aperture_sum_err" in photometry.colnames
                       else np.zeros(len(raw)))
            stats = ApertureStats(np.nan_to_num(data, nan=0.0), annulus,
                                  sigma_clip=SigmaClip(sigma=3.0, maxiters=5))
            sky = np.asarray(stats.median, dtype=float)
        except Exception as exc:
            self.logger.warning("Aperture photometry failed: %s", exc)
            return empty, empty, empty, empty

        flux = raw - np.nan_to_num(sky) * aperture.area
        flux_err = raw_err

        positive = np.isfinite(flux) & (flux > 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mag = np.where(positive, -2.5 * np.log10(np.abs(flux) / self.exposure), np.nan)
            mag_err = np.where(positive & (flux_err > 0),
                               _POGSON * flux_err / np.abs(flux), np.nan)
        if self.zero_point is not None:
            # The aperture magnitude gets the aperture zero point, not the
            # total-flux one: it is measured in that same aperture.
            mag = mag + (self.aperture_zero_point
                         if self.aperture_zero_point is not None else self.zero_point)
            mag_err = np.hypot(mag_err, self.zero_point_err)
        else:
            mag = np.full_like(mag, np.nan)
            mag_err = np.full_like(mag_err, np.nan)
        return flux, flux_err, mag, mag_err

    def _kron_magnitudes(self, cat):
        """Kron elliptical photometry -- SExtractor's ``AUTO``.

        The standard total magnitude for an extended source, and the thing the
        isophotal ``segment_flux`` is not: it adapts its aperture to each
        source's own light profile rather than stopping at a threshold.
        """
        empty = np.full(len(cat), np.nan)
        try:
            flux = np.asarray(cat.kron_flux, dtype=float)
            flux_err = np.asarray(cat.kron_fluxerr, dtype=float)
            radius = np.asarray(cat.kron_radius, dtype=float)
        except Exception as exc:
            self.logger.warning("Kron photometry unavailable: %s", exc)
            return empty, empty, empty, empty, empty

        positive = np.isfinite(flux) & (flux > 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mag = np.where(positive, -2.5 * np.log10(np.abs(flux) / self.exposure), np.nan)
            mag_err = np.where(positive & np.isfinite(flux_err) & (flux_err > 0),
                               _POGSON * flux_err / np.abs(flux), np.nan)
        if self.zero_point is not None:
            # Kron flux is already close to total, so it takes the total-flux
            # zero point rather than the aperture one.
            mag = mag + self.zero_point
            mag_err = np.hypot(mag_err, self.zero_point_err)
        else:
            mag = np.full_like(mag, np.nan)
            mag_err = np.full_like(mag_err, np.nan)
        return flux, flux_err, mag, mag_err, radius

    def _psf_magnitudes(self, data, err, tbl, fits_header):
        """PSF-fitted photometry, the matched filter for a point source.

        Uses an empirical PSF built from the frame's own isolated stars where
        one can be built, and an analytic Gaussian of the measured FWHM
        otherwise. A frame too sparse for either yields NaN, and classification
        then reports AMBIGUOUS rather than guessing.
        """
        empty = np.full(len(tbl), np.nan)
        if not self.config.phase3.steps.enabled("psf_photometry"):
            return empty, empty, empty
        try:
            from astropy.modeling.models import Gaussian2D
            from astropy.table import Table as _Table
            from photutils.psf import PSFPhotometry

            from cassa_photometry.psf import build_epsf, estimate_fwhm
        except ImportError:
            return empty, empty, empty

        try:
            model, info = build_epsf(data, estimate_fwhm(
                data, fwhm_guess=self.fwhm,
                threshold=self.config.phase4.detection_threshold,
                max_stars=self.config.phase4.max_stars,
                cutout=self.config.phase4.cutout,
            ))
            if model is None:
                self.logger.info(
                    "    No empirical PSF (%s); using an analytic Gaussian of the "
                    "measured FWHM.", info.get("reason"),
                )
                sigma = self.fwhm / 2.3548
                model = Gaussian2D(x_stddev=sigma, y_stddev=sigma)
                model.x_stddev.fixed = model.y_stddev.fixed = True

            init = _Table()
            init["x"] = np.asarray(tbl["xcentroid"], dtype=float)
            init["y"] = np.asarray(tbl["ycentroid"], dtype=float)
            fit_shape = int(2 * round(self.fwhm * 2) + 1)
            photometry = PSFPhotometry(
                model, fit_shape=max(fit_shape, 5),
                # Seeds each fit's flux from an aperture sum; without it
                # photutils refuses init_params that carry no flux column.
                aperture_radius=max(self.fwhm * self.config.phase3.aperture_r_factor, 3.0),
            )
            result = photometry(np.nan_to_num(data, nan=0.0), error=err, init_params=init)
        except Exception as exc:
            self.logger.warning("PSF photometry failed: %s", exc)
            return empty, empty, empty

        flux = np.asarray(result["flux_fit"], dtype=float)
        flux_err = np.asarray(result["flux_err"], dtype=float) \
            if "flux_err" in result.colnames else np.full(len(flux), np.nan)
        chi2 = (np.asarray(result["qfit"], dtype=float)
                if "qfit" in result.colnames else np.full(len(flux), np.nan))

        positive = np.isfinite(flux) & (flux > 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mag = np.where(positive, -2.5 * np.log10(np.abs(flux) / self.exposure), np.nan)
            mag_err = np.where(positive & np.isfinite(flux_err) & (flux_err > 0),
                               _POGSON * flux_err / np.abs(flux), np.nan)
        if self.zero_point is not None:
            mag = mag + self.zero_point
            mag_err = np.hypot(mag_err, self.zero_point_err)
        else:
            mag = np.full_like(mag, np.nan)
            mag_err = np.full_like(mag_err, np.nan)
        return mag, mag_err, chi2

    def _require_zp(self):
        if self.zero_point is None:
            raise RuntimeError("Zero point missing. Run calculate_local_zero_point first.")


def combine_zeropoints(zps, zp_errs, sigma=3.0):
    """Combine per-star zero points into a value with an uncertainty.

    Sigma-clips the offsets, then takes the inverse-variance weighted mean. The
    reported uncertainty is the larger of the formal weighted error and the
    standard error of the surviving scatter (a conservative choice).

    Returns
    -------
    (zero_point, zero_point_err, n_stars) : tuple(float, float, int)
    """
    zps = np.asarray(zps, dtype=float)
    zp_errs = np.asarray(zp_errs, dtype=float)
    keep = ~np.ma.getmaskarray(sigma_clip(zps, sigma=sigma, maxiters=5))
    zps, zp_errs = zps[keep], zp_errs[keep]

    weights = 1.0 / zp_errs ** 2
    zero_point = float(np.sum(weights * zps) / np.sum(weights))
    formal_err = np.sqrt(1.0 / np.sum(weights))
    scatter_err = np.std(zps, ddof=1) / np.sqrt(len(zps)) if len(zps) > 1 else formal_err
    return zero_point, float(max(formal_err, scatter_err)), int(len(zps))


def _ellipticity(tbl, cat):
    if "ellipticity" in tbl.colnames:
        return tbl["ellipticity"]
    if hasattr(cat, "ellipticity"):
        return cat.ellipticity
    if "eccentricity" in tbl.colnames:
        return tbl["eccentricity"]
    if hasattr(cat, "eccentricity"):
        return cat.eccentricity
    return np.full(len(tbl), np.nan)


def _segment_flags(segmap_data, labels, dq):
    """Bitwise-OR the DQ flags contained within each source segment."""
    if dq is None:
        return np.zeros(len(labels), dtype=np.int32)
    flags = ndimage.labeled_comprehension(
        dq, segmap_data, labels,
        lambda v: int(np.bitwise_or.reduce(v)) if len(v) else 0,
        np.int32, 0,
    )
    return np.asarray(flags, dtype=np.int32)


# --- Module helpers -----------------------------------------------------------

def _background_of(data, config):
    """One 2D background model, shared by the zero-point and catalog paths.

    The zero point used to detect against a single scalar median while the
    catalog used ``Background2D``: three different notions of "sky" in one
    module, only one of them spatially varying.
    """
    return Background2D(
        data, config.phase2.background_box,
        filter_size=config.phase2.background_filter,
        bkg_estimator=MedianBackground(),
    )


def _single_frame_exposure(header, logger=None):
    """The exposure an instrumental magnitude should be normalised by.

    A master is a weighted **mean**, so its pixels hold electrons per *single
    frame* -- not per ``TOT_EXP``. Dividing by ``TOT_EXP`` would be wrong by a
    factor of N. ``EXPMEAN`` is the true mean single-frame exposure, unrounded;
    ``EXPTIME`` is the fallback, and note that the group key quantises exposure
    to 5 s, so a 118 s frame recorded as 120 s carries a 2% (0.02 mag) error.
    """
    for card in ("EXPMEAN", "EXPTIME", "EXPOSURE"):
        try:
            value = float(header.get(card))
        except (TypeError, ValueError):
            continue
        if np.isfinite(value) and value > 0:
            if card != "EXPMEAN" and logger is not None:
                logger.info(
                    "    No EXPMEAN card; normalising by %s = %.4g s.", card, value
                )
            return value
    if logger is not None:
        logger.warning(
            "    No exposure time in the header; instrumental magnitudes are "
            "per-frame totals and the zero point is not comparable between epochs."
        )
    return 1.0


def _mutual_matches(det_coords, cat_coords, radius):
    """Mutually-nearest detection/catalog pairs within ``radius``.

    One-directional matching lets two detections claim the same catalog star --
    a deblended pair, or a star beside a galaxy -- and both then enter the zero
    point. Requiring the match to be nearest in *both* directions removes that,
    and costs one extra search.
    """
    if len(det_coords) == 0 or len(cat_coords) == 0:
        return []
    forward_idx, forward_d2d, _ = match_coordinates_sky(det_coords, cat_coords)
    reverse_idx, _, _ = match_coordinates_sky(cat_coords, det_coords)

    pairs = []
    for detection, (catalog_index, separation) in enumerate(
        zip(forward_idx, forward_d2d, strict=True)
    ):
        if separation > radius:
            continue
        if reverse_idx[catalog_index] != detection:
            continue  # some other detection is closer to this catalog star
        pairs.append((detection, int(catalog_index)))
    return pairs


def _aperture_has_bad_pixels(dq, aperture):
    """True for apertures containing a saturated, bad or cosmic-ray pixel."""
    from cassa_photometry.fits_utils import DQ_BAD_PIXEL, DQ_COSMIC_RAY, DQ_SATURATED

    bad = (np.asarray(dq).astype(np.int32)
           & (DQ_SATURATED | DQ_BAD_PIXEL | DQ_COSMIC_RAY)) != 0
    if not bad.any():
        return np.zeros(len(np.atleast_2d(aperture.positions)), dtype=bool)
    counts = aperture_photometry(bad.astype(float), aperture)["aperture_sum"]
    return np.asarray(counts, dtype=float) > 0


def _flux_radius(cat, fraction=0.5):
    """Half-light radius per source, a concentration measure."""
    try:
        return np.asarray(cat.fluxfrac_radius(fraction), dtype=float)
    except Exception:
        return np.full(len(cat), np.nan)


def _source_fwhm(tbl, cat):
    """Per-source FWHM, where photutils provides one."""
    for source in (tbl, cat):
        try:
            if hasattr(source, "colnames") and "fwhm" in source.colnames:
                return np.asarray(source["fwhm"], dtype=float)
            value = getattr(source, "fwhm", None)
            if value is not None:
                return np.asarray(value, dtype=float)
        except Exception:
            continue
    return np.full(len(tbl), np.nan)


def _has_flag(flags, which):
    """Sources carrying a given DQ flag."""
    from cassa_photometry.fits_utils import DQ_SATURATED

    bits = {"saturated": DQ_SATURATED}[which]
    return (np.asarray(flags, dtype=np.int64) & bits) != 0


def _near_edge(tbl, shape, fwhm, margin_fwhm=3.0):
    """Sources close enough to the border that their photometry is truncated."""
    margin = max(float(fwhm) * margin_fwhm, 5.0)
    x = np.asarray(tbl["xcentroid"], dtype=float)
    y = np.asarray(tbl["ycentroid"], dtype=float)
    ny, nx = shape
    return (x < margin) | (y < margin) | (x > nx - margin) | (y > ny - margin)
