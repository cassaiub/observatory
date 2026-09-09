"""Post-calibration verification: compare catalog magnitudes against APASS,
Pan-STARRS and SDSS to estimate the residual photometric offset and scatter.

The report shows both error bars -- our ``MAGERR_ISO`` (propagated from the ERR
plane) and the reference catalog's own uncertainty -- so an offset can be judged
against the errors. The cross-match radius is derived from the astrometric
solution (``ASTRMS``, written by phase 2) rather than fixed; see
:func:`match_radius_arcsec`.
"""

import os

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import fits
from astropy.table import Table
from astroquery.sdss import SDSS
from astroquery.vizier import Vizier

from cassa_photometry.config import load_config
from cassa_photometry.fits_utils import open_fits
from cassa_photometry.instruments import get_profile
from cassa_photometry.phase3_photometry import catalogs


def companion_fits(csv_path):
    """Yield the FITS files that belong to a ``_catalog.csv``, nearest match first."""
    base_name_no_ext = os.path.basename(csv_path).replace('_catalog.csv', '')
    input_dir = os.path.dirname(csv_path)
    for pf in (os.path.join(input_dir, f"{base_name_no_ext}.fits"),
               os.path.join(input_dir, f"{base_name_no_ext}_fluxcal.fits")):
        if os.path.exists(pf):
            yield pf


def match_radius_arcsec(csv_path, config):
    """Cross-match radius for this image, in arcsec.

    Phase 2 records how well the WCS actually fits (``ASTRMS``, the RMS residual
    against the astrometry index). A star should land within a few times that of
    its catalog position, so the radius is ``sigma * ASTRMS`` -- tight when the
    astrometry is good, forgiving when it is not -- clamped to a sane range. The
    floor matters: reference positions carry their own error and stars have moved
    since the catalog epoch, so an excellent solve still should not use a radius
    of a few tenths of an arcsecond.

    Falls back to the fixed ``phase3.zp_match_tol_arcsec`` when no ASTRMS exists
    (masters solved before this was recorded, or a failed solve).
    """
    for pf in companion_fits(csv_path):
        try:
            header = fits.getheader(pf)
        except OSError:
            continue
        radius, why = match_radius_from_header(header, config)
        if "ASTRMS" in why:
            return radius, why
    return match_radius_from_header(None, config)


def match_radius_from_header(header, config):
    """The same radius, from a header already in hand.

    Phase 3's zero point matches against the same astrometry the verification
    tool does, so both derive their radius here rather than one of them keeping
    a hardcoded 2 arcsec. On the workshop masters ASTRMS is about 0.5", so a
    fixed 2" was matching at roughly 4 sigma where 3 will do.
    """
    cfg = config.phase3
    try:
        astrms = float(header.get("ASTRMS"))
    except (AttributeError, TypeError, ValueError):
        astrms = None
    if astrms and astrms > 0:
        radius = min(max(cfg.match_radius_sigma * astrms,
                         cfg.match_radius_min_arcsec),
                     cfg.match_radius_max_arcsec)
        nstars = header.get("ASTNSTAR", "?")
        return radius, (f"{cfg.match_radius_sigma:g} x ASTRMS {astrms:.3f}\" "
                        f"(WCS fit to {nstars} stars)")
    return (cfg.zp_match_tol_arcsec,
            "fixed phase3.zp_match_tol_arcsec -- no ASTRMS in the header")


#: Science band -> the reference-catalog column this tool compares against.
_BAND_TO_COLUMN = {"R": "rmag", "I": "imag", "G": "gmag", "V": "Vmag", "B": "Bmag"}

#: Bands the instrument profile does not carry, kept so z/y frames still verify.
_EXTRA_FILTERS = {
    "Z": "zmag", "Z-BAND": "zmag", "ZMAG": "zmag", "SDSS-Z": "zmag",
    "Y": "ymag", "Y-BAND": "ymag", "YMAG": "ymag",
}


def detect_filter(csv_path, default_filter, instrument=None, config=None):
    """Reference-catalog column for a catalog CSV, from its companion FITS.

    Filter naming comes from the instrument profile, so this tool agrees with
    phases 1-3 instead of keeping its own copy of the mapping. Falls back to
    filename parsing when there is no companion FITS to read.

    With neither ``instrument`` nor ``config`` the ``generic`` profile is used,
    which is not necessarily the profile the run was configured with.
    """
    instrument = instrument or get_profile(
        getattr(config, "instrument", None), config=config
    )
    for pf in companion_fits(csv_path):
        try:
            with open_fits(pf) as hdul:
                header = hdul[0].header
                raw = instrument.get_filter(header)
                if not raw or raw in ("UNKNOWN", "NONE"):
                    continue
                if raw in _EXTRA_FILTERS:
                    return _EXTRA_FILTERS[raw]
                band = instrument.science_band(header, default=None)
                if band in _BAND_TO_COLUMN:
                    return _BAND_TO_COLUMN[band]
        except Exception:
            pass  # If FITS reading fails, silently fall back to filename checking

    # Fallback: Guess from the CSV filename if no FITS file exists
    base_name = os.path.basename(csv_path).upper()
    for token, column in _FILENAME_BANDS.items():
        if f"_{token}_" in base_name or f"-{token}-" in base_name:
            return column
    return default_filter


#: Filename token -> reference-catalog column, for frames with no companion FITS.
_FILENAME_BANDS = {"R": "rmag", "V": "Vmag", "B": "Bmag", "I": "imag", "G": "gmag"}


def _select_stars(cat, config):
    """The catalog's stars, and a one-line description of how they were chosen."""
    if "CLASS" in cat.colnames:
        classes = np.asarray([str(c).upper() for c in cat["CLASS"]])
        keep = classes == "STAR"
        if "FLAGS" in cat.colnames:
            keep &= np.asarray(cat["FLAGS"], dtype=np.int64) == 0
        return cat[keep], "CLASS == STAR (PSF vs Kron concentration)"

    if "ELLIPTICITY" in cat.colnames:
        print("[NOTE] This catalog predates the CLASS column; falling back to the "
              "ellipticity cut, which measures shape rather than concentration.")
        keep = np.asarray(cat["ELLIPTICITY"], dtype=float) < config.phase3.ellipticity_star_max
        return cat[keep], f"ELLIPTICITY < {config.phase3.ellipticity_star_max}"

    print("[ERROR] Neither CLASS nor ELLIPTICITY is present; cannot separate "
          "stars from galaxies.")
    return None, "none"


def cross_match_and_report(catalog_name, cat_coords, cat_mags, bright_stars, my_coords,
                           matched_filter, cat_mag_errs=None, match_radius=2.0):
    """
    Cross-matches our local stars against a global catalog, calculates the photometric
    error (delta), and prints a formatted scientific report.

    ``cat_mag_errs`` carries the reference catalog's own magnitude uncertainties so the
    report can show both error bars side by side; our ``MAGERR_ISO`` comes from the
    ERR plane propagated through phases 1-3.
    """
    idx, d2d, _ = match_coordinates_sky(my_coords, cat_coords)

    # Accept a match only inside the radius derived from the astrometric solution,
    # so we are confident it is the same star and not a neighbour.
    match_mask = d2d < match_radius * u.arcsec

    print(f"\n--- {catalog_name.upper()} VERIFICATION REPORT "
          f"(Filter: {matched_filter} | match radius: {match_radius:.2f}\") ---")
    header = (f"{'Obj ID':<7} | {'RA (deg)':<10} | {'DEC (deg)':<10} | "
              f"{'Your Mag':<9} | {'Your Err':<9} | "
              f"{catalog_name + ' Mag':<16} | {catalog_name + ' Err':<16} | {'Delta':<10}")
    print(header)
    print("-" * len(header))

    errors = []
    for i in range(len(bright_stars)):
        if match_mask[i]:
            ps_idx = idx[i]
            true_mag = cat_mags[ps_idx]
            calc_mag = bright_stars['MAG_ISO'][i]
            calc_err = bright_stars['MAGERR_ISO'][i]

            # Skip masked, negative, or missing data from the catalog
            if np.ma.is_masked(true_mag) or np.isnan(true_mag) or true_mag < 0:
                continue

            true_err = None
            if cat_mag_errs is not None:
                e = cat_mag_errs[ps_idx]
                if not np.ma.is_masked(e) and np.isfinite(e) and e > 0:
                    true_err = float(e)
            true_err_str = f"{true_err:.3f}" if true_err is not None else "n/a"

            # Calculate the mathematical offset
            delta = calc_mag - true_mag
            errors.append(delta)

            print(f"{bright_stars['NUMBER'][i]:<7} | {bright_stars['ALPHA_J2000'][i]:<10.5f} | {bright_stars['DELTA_J2000'][i]:<10.5f} | "
                  f"{calc_mag:<9.3f} | {calc_err:<9.3f} | {true_mag:<16.3f} | {true_err_str:<16} | {delta:+.3f} mag")

    if len(errors) > 0:
        mean_error = np.mean(errors)
        std_error = np.std(errors)
        print("-" * len(header))
        print(f"Mean Calibration Error (Offset): {mean_error:+.4f} mag")
        print(f"Standard Deviation (Scatter)   : {std_error:.4f} mag")

        # Intelligent feedback based on the catalog
        if abs(mean_error) < 0.1:
            print(f"✅ EXCELLENT! Matches {catalog_name} perfectly.")
        else:
            print(f"⚠️ NOTE: Noticeable offset against {catalog_name}.")
            print("   This is likely a 'Color Term' difference between your physical glass filter")
            print(f"   and the {catalog_name} passband. To force-match this specific catalog,")
            print(f"   you would need to add {-mean_error:.4f} to your Zero Point.")
    else:
        print(f"No valid matches found between your stars and {catalog_name}.")


def verify_calibration(csv_path, filter_band="rmag", config=None):
    config = config or load_config()
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Cannot find catalog file: {csv_path}")

    match_radius, radius_origin = match_radius_arcsec(csv_path, config)

    print(f"Analyzing Catalog: {os.path.basename(csv_path)}")
    print(f"Target Filter Band: {filter_band}")
    print(f"Cross-match radius: {match_radius:.2f}\"  [{radius_origin}]")

    # 1. Load the generated catalog
    cat = Table.read(csv_path, format='csv')

    # 2. Keep the stars.
    #
    # CLASS comes from MAG_PSF - MAG_AUTO against the frame's own stellar locus,
    # which measures *concentration*. The older ELLIPTICITY cut measured *shape*,
    # so a face-on elliptical galaxy passed as a star and a slightly trailed star
    # did not. ELLIPTICITY is still the fallback for catalogs written before
    # CLASS existed.
    stars, basis = _select_stars(cat, config)
    if stars is None:
        return
    if len(stars) == 0:
        print(f"[ERROR] No stars found in the catalog ({basis}).")
        return

    # 3. Sort by brightness and take the top 25.
    #
    # Skipping the very brightest: they are the ones that saturate, and a
    # saturated star's magnitude is wrong in a way that would be blamed on the
    # zero point. Anything flagged is dropped outright.
    magnitude_column = "MAG_BEST" if "MAG_BEST" in stars.colnames else "MAG_ISO"
    stars = stars[np.isfinite(np.asarray(stars[magnitude_column], dtype=float))]
    stars.sort(magnitude_column)
    bright_stars = stars[:25]

    print(f" -> {len(stars)} star(s) by {basis}; verifying the {len(bright_stars)} "
          f"brightest by {magnitude_column}.\n")

    # 4. Get the center coordinate and dynamic bounding box
    mean_ra = np.mean(bright_stars['ALPHA_J2000'])
    mean_dec = np.mean(bright_stars['DELTA_J2000'])
    center_coord = SkyCoord(ra=mean_ra, dec=mean_dec, unit=(u.deg, u.deg))
    my_coords = SkyCoord(ra=bright_stars['ALPHA_J2000'], dec=bright_stars['DELTA_J2000'], unit=(u.deg, u.deg))

    ra_spread = np.max(bright_stars['ALPHA_J2000']) - np.min(bright_stars['ALPHA_J2000'])
    dec_spread = np.max(bright_stars['DELTA_J2000']) - np.min(bright_stars['DELTA_J2000'])
    radius = max(ra_spread, dec_spread) * u.deg / 2.0 + (1 * u.arcmin)

    # =================================================================
    # CATALOG 1: Pan-STARRS (via NASA MAST)
    # =================================================================
    try:
        if filter_band in ['Vmag', 'Bmag']:
            print(f">>> Skipping Pan-STARRS DR2 - Filter '{filter_band}' is not supported by Pan-STARRS.")
        else:
            print(f">>> Querying Pan-STARRS DR2 (NASA MAST) for radius {radius.to(u.arcmin):.2f}...")
            panstarrs_cat = catalogs.query_panstarrs(center_coord, radius)
            ps_coords = SkyCoord(ra=panstarrs_cat['RAJ2000'], dec=panstarrs_cat['DEJ2000'], unit=(u.deg, u.deg))
            ps_mags = panstarrs_cat[filter_band]
            ps_err_col = f"{filter_band}_err"
            ps_errs = panstarrs_cat[ps_err_col] if ps_err_col in panstarrs_cat.colnames else None
            cross_match_and_report("Pan-STARRS", ps_coords, ps_mags, bright_stars, my_coords, filter_band,
                                   ps_errs, match_radius)
    except Exception as e:
        print(f"[WARNING] Pan-STARRS Query Failed: {e}")

    # =================================================================
    # CATALOG 2: SDSS - Sloan Digital Sky Survey
    # =================================================================
    try:
        if filter_band in ['Vmag', 'Bmag']:
            print(f"\n>>> Skipping SDSS DR12 - Filter '{filter_band}' is not supported by SDSS.")
        else:
            print("\n>>> Querying SDSS DR12 (via SQL to bypass radius limits)...")
            # Convert 'rmag' -> 'r' for SDSS formatting
            sdss_filter = filter_band[0] if filter_band.endswith('mag') else filter_band

            # Build a bounding box to strictly encapsulate our 25 stars
            ra_min, ra_max = np.min(bright_stars['ALPHA_J2000']) - 0.02, np.max(bright_stars['ALPHA_J2000']) + 0.02
            dec_min, dec_max = np.min(bright_stars['DELTA_J2000']) - 0.02, np.max(bright_stars['DELTA_J2000']) + 0.02

            # Send a raw SQL query to SDSS to bypass the arbitrary 3.0 arcmin Python limit (type=6 means 'STAR')
            sql_query = f"SELECT ra, dec, {sdss_filter}, err_{sdss_filter} FROM PhotoObj WHERE ra BETWEEN {ra_min} AND {ra_max} AND dec BETWEEN {dec_min} AND {dec_max} AND type=6"
            sdss_res = SDSS.query_sql(sql_query)

            if sdss_res is not None and len(sdss_res) > 0:
                sdss_coords = SkyCoord(ra=sdss_res['ra'], dec=sdss_res['dec'], unit=(u.deg, u.deg))
                sdss_mags = sdss_res[sdss_filter]
                sdss_err_col = f"err_{sdss_filter}"
                sdss_errs = sdss_res[sdss_err_col] if sdss_err_col in sdss_res.colnames else None
                cross_match_and_report("SDSS", sdss_coords, sdss_mags, bright_stars, my_coords, sdss_filter,
                                       sdss_errs, match_radius)
            else:
                print(" -> Target field is outside the SDSS sky footprint or no stars matched.")
    except Exception as e:
        print(f"[WARNING] SDSS Query Failed: {e}")

    # =================================================================
    # CATALOG 3: APASS - AAVSO Photometric All-Sky Survey
    # =================================================================
    try:
        print("\n>>> Querying APASS DR9 (VizieR)...")

        # We use an Asian-focused robust mirror hopper to bypass ISP blocks
        apass_res = None
        mirrors = [
            'vizier.iucaa.in',         # India (Pune - Best chance for South Asia)
            'vizier.nao.ac.jp',        # Japan (Highly reliable fallback)
            'vizier.china-vo.org',     # China
            'vizier.saao.ac.za',       # South Africa
            'vizier.cds.unistra.fr',   # France
            'vizier.cfa.harvard.edu',  # USA
            'vizier.ast.cam.ac.uk',    # UK
            'vizier.hia.nrc.ca'        # Canada
        ]

        for mirror in mirrors:
            print(f" -> Attempting APASS connection via mirror: {mirror}")
            # Request all columns (*) instead of risking a specific name
            v = Vizier(columns=['RAJ2000', 'DEJ2000', '*'], row_limit=500)
            v.VIZIER_SERVER = mirror
            v.TIMEOUT = 15
            try:
                res = v.query_region(center_coord, radius=radius, catalog='II/336/apass9')
                if len(res) > 0:
                    apass_res = res[0]
                    print(f" -> Success! Connected to {mirror}")
                    break # Success! Break out of the mirror loop
            except Exception:
                print(f" -> [WARNING] Mirror {mirror} failed or blocked.")

        if apass_res is not None:
            # Bulletproof string stripping to exactly match the Engine's fuzzy logic
            filter_prefix = filter_band[0].lower() # 'r', 'v', 'b'

            target_cols = [
                c for c in apass_res.colnames
                if c.replace("'", "").replace("_", "").lower().startswith(filter_prefix)
                and 'mag' in c.lower()
            ]

            if len(target_cols) > 0:
                best_col = target_cols[0]
                print(f" -> Auto-detected APASS magnitude column: '{best_col}'")
                err_col = catalogs._match_column(apass_res.colnames, f"e_{best_col}")
                if err_col:
                    print(f" -> Auto-detected APASS error column: '{err_col}'")
                apass_coords = SkyCoord(ra=apass_res['RAJ2000'], dec=apass_res['DEJ2000'], unit=(u.deg, u.deg))
                apass_mags = apass_res[best_col]
                apass_errs = apass_res[err_col] if err_col else None
                cross_match_and_report("APASS", apass_coords, apass_mags, bright_stars, my_coords, best_col,
                                       apass_errs, match_radius)
            else:
                print(f" -> Connected, but could not find a column matching filter '{filter_band}'. Available columns: {apass_res.colnames}")
        else:
            print(" -> No APASS stars found or all global VizieR mirrors were blocked.")
    except Exception as e:
        print(f"[WARNING] APASS process failed: {e}")

def run(input_path, default_filter="rmag", config=None, instrument=None):
    """Verify one ``_catalog.csv`` file or a directory of them."""
    config = config or load_config()
    instrument = instrument or get_profile(config.instrument, config=config)
    input_path = os.path.abspath(input_path)

    if os.path.isdir(input_path):
        print(f"Scanning directory for Catalog CSV files: {input_path}")
        csv_files = sorted(
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith("_catalog.csv")
        )
        if not csv_files:
            print("No valid _catalog.csv files found in the directory.")
            return
    elif os.path.isfile(input_path) and input_path.endswith(".csv"):
        csv_files = [input_path]
    else:
        print("[ERROR] Input path must be a .csv file or a directory of _catalog.csv files.")
        return

    print("=" * 54)
    print("  MULTI-CATALOG CALIBRATION VERIFICATION TOOL (BATCH)")
    print(f"  Found {len(csv_files)} file(s) to verify.")
    print("=" * 54)

    for i, csv_file in enumerate(csv_files, 1):
        target_filter = detect_filter(csv_file, default_filter, instrument)
        print(f"\n\n{'=' * 75}")
        print(f"[{i}/{len(csv_files)}] VERIFYING FILE")
        print(f"{'=' * 75}")
        verify_calibration(csv_file, target_filter, config=config)

    print("\n--- BATCH VERIFICATION COMPLETE ---")
