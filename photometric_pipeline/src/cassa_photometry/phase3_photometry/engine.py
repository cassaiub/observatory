"""Photometry engine: zero point, flux calibration, and error-carrying catalogs.

Every measurement now carries an uncertainty. The zero point is a sigma-clipped,
inverse-variance weighted combination of per-star offsets and comes with a
``MAGZERR``; the output catalog exposes ``Flux_Error``, ``Mag_Error``, ``SNR``
and DQ ``FLAGS``. The ``ERR`` plane written by phases 1-2 drives the photometric
errors; if it is missing (legacy input) it is reconstructed from the background.
"""

import numpy as np
import astropy.units as u
from scipy import ndimage
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.convolution import Gaussian2DKernel, convolve

from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
from photutils.utils import calc_total_error

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.fits_utils import read_mef, write_mef
from cassa_photometry.phase3_photometry.catalogs import fetch_reference_catalog

# 2.5 / ln(10): converts a fractional flux error into a magnitude error.
_POGSON = 2.5 / np.log(10.0)


class UniversalPhotometryEngine:
    def __init__(self, fwhm_estimate=None, detection_threshold=None, config=None, logger=None):
        self.config = config or load_config()
        cfg = self.config.phase3
        self.fwhm = fwhm_estimate if fwhm_estimate is not None else cfg.fwhm
        self.threshold = detection_threshold if detection_threshold is not None else cfg.detection_threshold
        self.logger = logger or get_logger("cassa_photometry")
        self.zero_point = None
        self.zero_point_err = None
        self.n_zp_stars = 0

    # -- helpers --------------------------------------------------------------
    def _load(self, path):
        """Return ``(sci, err, dq, header)``; reconstruct ERR if it is missing."""
        sci, err, dq, header = read_mef(path)
        if err is None:
            _, _, std = sigma_clipped_stats(sci, sigma=3.0)
            err = calc_total_error(np.clip(sci, 0, None), np.full_like(sci, std), effective_gain=1.0)
            self.logger.warning("No ERR plane found; reconstructed background+Poisson errors.")
        return sci, err, dq, header

    # -- zero point -----------------------------------------------------------
    def calculate_local_zero_point(self, master_science_fits, science_band="R"):
        """Compute a filter-wise zero point with uncertainty (differential photometry)."""
        cfg = self.config.phase3
        self.logger.info("Calculating local zero point (differential photometry)...")
        data, err, _, header = self._load(master_science_fits)
        wcs = WCS(header)

        ny, nx = data.shape
        center = wcs.pixel_to_world(nx / 2, ny / 2)
        pixel_scale = np.abs(wcs.pixel_scale_matrix[0, 0]) * u.deg
        search_radius = (nx / 2) * pixel_scale

        catalog = fetch_reference_catalog(center, search_radius, science_band, self.logger)
        cat_coords = SkyCoord(ra=catalog["RAJ2000"], dec=catalog["DEJ2000"], unit=(u.deg, u.deg))

        _, median, std = sigma_clipped_stats(data, sigma=3.0)
        finder = DAOStarFinder(fwhm=self.fwhm, threshold=self.threshold * std)
        sources = finder.find_stars(data - median)
        if sources is None or len(sources) == 0:
            raise ValueError("No stars detected to determine a zero point.")

        det_coords = wcs.pixel_to_world(sources["xcentroid"], sources["ycentroid"])
        idx, d2d, _ = match_coordinates_sky(det_coords, cat_coords)
        matched = np.where(d2d < cfg.zp_match_tol_arcsec * u.arcsec)[0]

        zps, zp_errs = [], []
        for i in matched:
            cat_mag = float(catalog["Ref_Mag"][idx[i]])
            cat_mag_err = float(catalog["Ref_Mag_Err"][idx[i]])
            pos = (sources["xcentroid"][i], sources["ycentroid"][i])
            aperture = CircularAperture(pos, r=self.fwhm * cfg.aperture_r_factor)
            annulus = CircularAnnulus(pos, r_in=self.fwhm * cfg.annulus_in_factor,
                                      r_out=self.fwhm * cfg.annulus_out_factor)

            phot = aperture_photometry(data, [aperture, annulus], error=err)
            bkg_mean = phot["aperture_sum_1"][0] / annulus.area
            flux = phot["aperture_sum_0"][0] - bkg_mean * aperture.area
            # Flux error: source aperture error + subtracted-background error, in quadrature.
            flux_err = np.sqrt(
                phot["aperture_sum_err_0"][0] ** 2
                + (aperture.area / annulus.area) ** 2 * phot["aperture_sum_err_1"][0] ** 2
            )
            if flux <= 0 or not np.isfinite(flux_err) or flux_err <= 0:
                continue
            inst_mag = -2.5 * np.log10(flux)
            inst_mag_err = _POGSON * flux_err / flux
            zps.append(cat_mag - inst_mag)
            zp_errs.append(np.hypot(inst_mag_err, cat_mag_err))

        if len(zps) == 0:
            raise ValueError("Could not match any stars to the catalog for a zero point.")

        self.zero_point, self.zero_point_err, self.n_zp_stars = combine_zeropoints(
            np.array(zps), np.array(zp_errs), sigma=cfg.zp_sigma_clip,
        )
        self.logger.info(f"Zero point ({science_band}-band): {self.zero_point:.4f} "
                         f"+/- {self.zero_point_err:.4f} mag from {self.n_zp_stars} stars.")
        return self.zero_point

    # -- flux-calibrated image ------------------------------------------------
    def export_flux_calibrated_image(self, target_fits, output_filename):
        """Write a SCI/ERR image in Janskys (propagating the error plane)."""
        self._require_zp()
        conversion = self.config.phase3.ab_flux_zero_jy * 10 ** (-self.zero_point / 2.5)
        self.logger.info(f"Writing flux-calibrated image: {output_filename}")

        data, err, dq, header = self._load(target_fits)
        header["BUNIT"] = "Jy"
        header["MAGZERO"] = (self.zero_point, "Photometric zero point [mag]")
        header["MAGZERR"] = (self.zero_point_err, "Zero point uncertainty [mag]")
        header["NZPSTARS"] = (self.n_zp_stars, "Number of stars used for zero point")
        header["FLUXCAL"] = (conversion, "Jy per data unit")
        write_mef(output_filename, sci=data * conversion,
                  err=(None if err is None else err * conversion), dq=dq, header=header,
                  history=f"Flux calibrated to Jy using ZP {self.zero_point:.4f}")

    # -- source catalog -------------------------------------------------------
    def generate_full_catalog(self, target_fits, output_csv, output_segmap=None):
        """Detect sources and write a catalog with flux/magnitude uncertainties."""
        self._require_zp()
        cfg = self.config.phase3
        box = self.config.phase2.background_box
        self.logger.info("Generating full object catalog and segmentation map...")
        data, err, dq, header = self._load(target_fits)
        wcs = WCS(header)

        bkg = Background2D(data, box, filter_size=(3, 3), bkg_estimator=MedianBackground())
        data_sub = data - bkg.background
        threshold = self.threshold * bkg.background_rms

        kernel = Gaussian2DKernel(x_stddev=self.fwhm / 2.35)
        kernel.normalize()
        convolved = convolve(data_sub, kernel)

        segmap = detect_sources(convolved, threshold, npixels=cfg.detect_npixels)
        if segmap is None:
            self.logger.warning("No sources detected.")
            return
        segmap = deblend_sources(data_sub, segmap, npixels=cfg.detect_npixels,
                                 nlevels=cfg.deblend_nlevels, contrast=cfg.deblend_contrast)
        self.logger.info(f"Detected {segmap.nlabels} objects after deblending.")

        if output_segmap:
            seg_header = header.copy()
            seg_header["BUNIT"] = "ID"
            write_mef(output_segmap, sci=segmap.data.astype(np.int32), header=seg_header,
                      history="Segmentation map from UniversalPhotometryEngine")

        cat = SourceCatalog(data_sub, segmap, error=err, wcs=wcs)
        tbl = cat.to_table()

        flux = np.asarray(tbl["segment_flux"], dtype=float)
        flux_err = np.asarray(tbl["segment_fluxerr"], dtype=float)
        valid = flux > 0
        tbl, cat = tbl[valid], cat[valid]
        flux, flux_err = flux[valid], flux_err[valid]

        mag = -2.5 * np.log10(flux) + self.zero_point
        mag_err = np.sqrt((_POGSON * flux_err / flux) ** 2 + self.zero_point_err ** 2)

        out = Table()
        out["ID"] = tbl["label"]
        out["RA_deg"] = cat.sky_centroid.ra.deg
        out["Dec_deg"] = cat.sky_centroid.dec.deg
        out["X_pix"] = tbl["xcentroid"]
        out["Y_pix"] = tbl["ycentroid"]
        out["Instrumental_Flux"] = flux
        out["Flux_Error"] = flux_err
        out["Absolute_Mag"] = mag
        out["Mag_Error"] = mag_err
        out["SNR"] = np.where(flux_err > 0, flux / flux_err, np.nan)
        out["Area_pixels"] = tbl["area"]
        out["Ellipticity"] = _ellipticity(tbl, cat)
        out["FLAGS"] = _segment_flags(segmap.data, np.asarray(tbl["label"]), dq)

        out.write(output_csv, format="csv", overwrite=True)
        self.logger.info(f"Catalog saved to {output_csv} ({len(out)} sources).")

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
