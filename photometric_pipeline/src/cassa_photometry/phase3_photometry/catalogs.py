"""Reference-catalog access for photometric zero-point calibration.

Queries APASS (VizieR), Pan-STARRS (MAST) and SDSS with a filter-aware priority
cascade, and -- new for the error budget -- also returns the catalog magnitude
*uncertainty* (``Ref_Mag_Err``) so it can be folded into the zero-point error.
"""

import urllib.request

import numpy as np
import astropy.units as u
from astropy.table import Table
from astroquery.vizier import Vizier
from astroquery.sdss import SDSS

from cassa_photometry.logging_utils import get_logger

_APASS_MIRRORS = [
    "vizier.cfa.harvard.edu",
    "vizier.cds.unistra.fr",
    "vizier.nao.ac.jp",
    "vizier.iucaa.in",
    "vizier.china-vo.org",
]

# science band -> per-catalog (magnitude column, error column) hints.
_BAND_COLUMNS = {
    "R": {"APASS": ("r'mag", "e_r'mag"), "Pan-STARRS": ("rmag", "rmag_err"), "SDSS": ("r", "err_r")},
    "G": {"APASS": ("g'mag", "e_g'mag"), "Pan-STARRS": ("gmag", "gmag_err"), "SDSS": ("g", "err_g")},
    "I": {"APASS": ("i'mag", "e_i'mag"), "Pan-STARRS": ("imag", "imag_err"), "SDSS": ("i", "err_i")},
    "V": {"APASS": ("Vmag", "e_Vmag")},
    "B": {"APASS": ("Bmag", "e_Bmag")},
}

# Default magnitude uncertainty when a catalog omits an error column.
_DEFAULT_MAG_ERR = 0.03


def query_apass(center_coord, radius):
    """Query APASS DR9 via VizieR, hopping mirrors to bypass regional blocks."""
    for mirror in _APASS_MIRRORS:
        v = Vizier(columns=["RAJ2000", "DEJ2000", "*"], row_limit=1000)
        v.TIMEOUT = 15
        v.VIZIER_SERVER = mirror
        try:
            result = v.query_region(center_coord, radius=radius, catalog="II/336/apass9")
            if len(result) > 0:
                return result[0]
        except Exception:
            continue
    return None


def query_panstarrs(center_coord, radius):
    """Query Pan-STARRS DR2 mean object catalog via the MAST REST API."""
    try:
        ra_deg, dec_deg = center_coord.ra.deg, center_coord.dec.deg
        radius_deg = radius.to(u.deg).value
        url = (f"https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
               f"?ra={ra_deg}&dec={dec_deg}&radius={radius_deg}")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as response:
            csv_text = response.read().decode("utf-8")

        lines = [line.strip() for line in csv_text.strip().split("\n") if line.strip()]
        if len(lines) < 2:
            return None
        header = [h.strip() for h in lines[0].split(",")]

        def col(name):
            return header.index(name) if name in header else -1

        idx = {
            "ra": col("raMean"), "dec": col("decMean"),
            "gmag": col("gMeanPSFMag"), "rmag": col("rMeanPSFMag"), "imag": col("iMeanPSFMag"),
            "gmag_err": col("gMeanPSFMagErr"), "rmag_err": col("rMeanPSFMagErr"),
            "imag_err": col("iMeanPSFMagErr"),
        }
        if idx["ra"] < 0 or idx["dec"] < 0 or idx["rmag"] < 0:
            return None

        rows = {k: [] for k in idx}
        for line in lines[1:]:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) <= max(idx["ra"], idx["dec"], idx["rmag"]):
                continue
            try:
                ra_v, dec_v = float(parts[idx["ra"]]), float(parts[idx["dec"]])
            except ValueError:
                continue
            rows["ra"].append(ra_v)
            rows["dec"].append(dec_v)
            for key in ("gmag", "rmag", "imag", "gmag_err", "rmag_err", "imag_err"):
                rows[key].append(_parse_float(parts, idx[key]))

        return Table(
            [rows["ra"], rows["dec"], rows["gmag"], rows["rmag"], rows["imag"],
             rows["gmag_err"], rows["rmag_err"], rows["imag_err"]],
            names=("RAJ2000", "DEJ2000", "gmag", "rmag", "imag",
                   "gmag_err", "rmag_err", "imag_err"),
        )
    except Exception:
        return None


def query_sdss(center_coord, radius):
    """Query SDSS DR12 point sources via raw SQL (bypasses the radius cap)."""
    try:
        r = radius.to(u.deg).value
        ra0, dec0 = center_coord.ra.deg, center_coord.dec.deg
        sql = (
            "SELECT ra as RAJ2000, dec as DEJ2000, u, g, r, i, z, "
            "err_u, err_g, err_r, err_i, err_z FROM PhotoObj "
            f"WHERE ra BETWEEN {ra0 - r} AND {ra0 + r} "
            f"AND dec BETWEEN {dec0 - r} AND {dec0 + r} AND type=6"
        )
        return SDSS.query_sql(sql)
    except Exception:
        return None


def fetch_reference_catalog(center_coord, radius, science_band, logger=None):
    """Return a standardised reference catalog for ``science_band``.

    Returns
    -------
    astropy.table.Table
        Columns ``RAJ2000``, ``DEJ2000``, ``Ref_Mag``, ``Ref_Mag_Err``.

    Raises
    ------
    ValueError
        If ``science_band`` is unsupported.
    RuntimeError
        If every prioritised catalog fails.
    """
    logger = logger or get_logger("cassa_photometry")
    band = science_band.upper().strip()
    if band not in _BAND_COLUMNS:
        raise ValueError(f"Unsupported science band: {science_band!r}")

    queries = {"APASS": query_apass, "Pan-STARRS": query_panstarrs, "SDSS": query_sdss}
    priority = [(name, cols[0], cols[1], queries[name])
                for name, cols in _BAND_COLUMNS[band].items()]

    logger.info(f"Catalog priority for band '{band}': " + ", ".join(p[0] for p in priority))

    for cat_name, mag_hint, err_hint, query_func in priority:
        logger.info(f"[ATTEMPTING] {cat_name}...")
        raw = query_func(center_coord, radius)
        if raw is None or len(raw) == 0:
            logger.warning(f"[FAILED] {cat_name} returned no stars.")
            continue

        mag_col = _match_column(raw.colnames, mag_hint, require="mag")
        if not mag_col:
            logger.warning(f"[FAILED] {cat_name}: no column matching '{mag_hint}'.")
            continue
        err_col = _match_column(raw.colnames, err_hint)

        valid = ~np.isnan(np.asarray(raw[mag_col], dtype=float))
        if hasattr(raw[mag_col], "mask"):
            valid &= ~raw[mag_col].mask
        clean = raw[valid]
        if len(clean) == 0:
            logger.warning(f"[FAILED] {cat_name}: no valid {mag_hint} values.")
            continue

        std = Table()
        std["RAJ2000"] = clean["RAJ2000"]
        std["DEJ2000"] = clean["DEJ2000"]
        std["Ref_Mag"] = np.asarray(clean[mag_col], dtype=float)
        if err_col:
            errs = np.asarray(clean[err_col], dtype=float)
            errs[~np.isfinite(errs) | (errs <= 0)] = _DEFAULT_MAG_ERR
            std["Ref_Mag_Err"] = errs
        else:
            std["Ref_Mag_Err"] = np.full(len(clean), _DEFAULT_MAG_ERR)
        logger.info(f"[SUCCESS] {len(std)} reference stars from {cat_name} "
                    f"(mag='{mag_col}', err='{err_col or _DEFAULT_MAG_ERR}').")
        return std

    raise RuntimeError(f"All prioritised catalogs failed for band '{band}'.")


def _parse_float(parts, idx):
    if idx == -1 or idx >= len(parts):
        return np.nan
    try:
        val = float(parts[idx])
        return val if val > -99.0 else np.nan
    except ValueError:
        return np.nan


def _match_column(colnames, hint, require=None):
    """Fuzzy-match a column name, ignoring apostrophes/underscores/case."""
    if not hint:
        return None
    if hint in colnames:
        return hint
    target = hint.replace("'", "").replace("_", "").lower()
    for c in colnames:
        cclean = c.replace("'", "").replace("_", "").lower()
        if target in cclean and (require is None or require in c.lower()):
            return c
    return None
