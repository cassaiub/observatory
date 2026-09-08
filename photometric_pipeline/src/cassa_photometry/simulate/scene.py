"""The sky half of the forward model: a star field and a galaxy, in e-/s.

Every choice here exists to make some pipeline defect *detectable*:

* **A Moffat PSF, not a Gaussian.** Enclosed flux inside 2xFWHM is 99.998% for a
  Gaussian and 93.4% for a Moffat with beta = 2.5. A Gaussian scene therefore has
  an aperture correction of essentially zero, which is exactly why a synthetic
  test suite can pass while the real aperture correction is missing.
* **A PSF that broadens toward the field edge.** An 8-inch f/5 Newtonian has
  coma, so a *scalar* aperture correction is right at the centre and wrong in the
  corners. A constant-PSF scene cannot show that.
* **A range of magnitudes crossing saturation.** The brightest catalogue stars
  are the ones that saturate and the ones inverse-variance weighting trusts most.
* **Colours that differ star to star**, so a colour term is measurable rather
  than absorbed into the zero point.
"""

import os

import numpy as np
from astropy.wcs import WCS

#: Star fields shipped with the package, as real Gaia DR3 positions and
#: synthetic Johnson-Cousins magnitudes. Using a real field is what lets
#: ``solve-field`` actually solve a simulated frame and lets the reference-catalog
#: cross-match find the stars -- so an end-to-end test exercises astrometry and
#: photometric calibration rather than only the arithmetic in between.
FIELD_DIR = os.path.join(os.path.dirname(__file__), "data")

#: Catalog column holding each band's magnitude, in the Johnson-Cousins system
#: the ``B``/``V``/``R``/``I`` filters are actually on.
CATALOG_COLUMNS = {"B": "b_jkc_mag", "V": "v_jkc_mag",
                   "R": "r_jkc_mag", "I": "i_jkc_mag"}

#: Instrumental zero point per band: the magnitude of a source giving 1 e-/s.
#: Representative of a 203 mm aperture at ~40% total throughput, and the value
#: the pipeline should recover once exposure normalisation is correct.
ZERO_POINTS = {"B": 19.2, "V": 20.0, "R": 20.1, "I": 19.4}

#: First-order atmospheric extinction, magnitudes per airmass.
EXTINCTION = {"B": 0.25, "V": 0.15, "R": 0.10, "I": 0.06}

#: Sky surface brightness per band, mag/arcsec^2 at zenith. Suburban.
SKY_BRIGHTNESS = {"B": 19.5, "V": 18.8, "R": 18.3, "I": 17.8}

#: Moffat beta. 2.5 is the usual atmospheric value and gives realistic wings.
MOFFAT_BETA = 2.5


def moffat_alpha(fwhm, beta=MOFFAT_BETA):
    """Moffat core width from a FWHM."""
    return fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))


def enclosed_fraction(radius, fwhm, beta=MOFFAT_BETA):
    """Fraction of a Moffat's total flux inside ``radius`` pixels.

    Used by the tests to state the expected aperture correction exactly, rather
    than measuring it from the same code that produced the image.
    """
    alpha = moffat_alpha(fwhm, beta)
    return 1.0 - (1.0 + (radius / alpha) ** 2) ** (1.0 - beta)


class Scene:
    """A patch of sky: sources with known magnitudes, plus a galaxy.

    The truth catalogue this exposes is the point of the whole module -- every
    photometric assertion downstream compares against it.
    """

    def __init__(self, shape, pixel_scale, centre_ra, centre_dec, rng=None,
                 n_stars=250, galaxy=True, coma_strength=0.6, field=None):
        self.shape = tuple(shape)
        self.pixel_scale = float(pixel_scale)
        self.rng = rng if rng is not None else np.random.default_rng(0)
        self.coma_strength = float(coma_strength)
        self.wcs = self._build_wcs(centre_ra, centre_dec)
        self.field = field
        self.stars = (
            self._build_stars_from_field(field, n_stars) if field
            else self._build_stars_synthetic(n_stars)
        )
        self.galaxy = self._build_galaxy() if galaxy else None

    # --- Geometry -------------------------------------------------------------
    def _build_wcs(self, centre_ra, centre_dec):
        ny, nx = self.shape
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [nx / 2.0 + 0.5, ny / 2.0 + 0.5]
        wcs.wcs.crval = [centre_ra, centre_dec]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        scale_deg = self.pixel_scale / 3600.0
        wcs.wcs.cdelt = [-scale_deg, scale_deg]
        return wcs

    def _relative_field_radius(self, x, y):
        """Distance from the field centre, 0 at centre and 1 at the corner."""
        ny, nx = self.shape
        cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
        return np.hypot(x - cx, y - cy) / np.hypot(cx, cy)

    def local_fwhm(self, x, y, base_fwhm):
        """Delivered FWHM at a position: broader off-axis, because of coma."""
        return base_fwhm * (1.0 + self.coma_strength * self._relative_field_radius(x, y) ** 2)

    # --- Sources --------------------------------------------------------------
    def _build_stars_from_field(self, field, n_stars):
        """Real stars from a shipped Gaia field, plus a synthetic faint tail.

        Only sources landing on the detector are kept. Magnitudes are Gaia DR3
        synthetic photometry in Johnson-Cousins, so the simulated frames are on
        the same system the filters claim to be -- which is what makes a
        recovered zero point comparable to truth without a colour term standing
        in the way.

        **Why a synthetic tail is added on top.** Gaia synthetic photometry
        requires an XP spectrum, so the catalogue stops around V = 18. A frame
        built from it alone is therefore *sparser than the real sky*, which makes
        alignment, detection and crowding all easier than they are in practice --
        the failure modes worth testing would be simulated away. The faint tail
        continues the observed counts below the catalogue limit. Those stars are
        equally truth-known; they simply have no Gaia source id.
        """
        import pandas as pd

        path = field if os.path.sep in str(field) else os.path.join(
            FIELD_DIR, f"field_{field}.csv"
        )
        table = pd.read_csv(path)

        x, y = self.wcs.all_world2pix(table["ra"].to_numpy(), table["dec"].to_numpy(), 0)
        ny, nx = self.shape
        on_chip = (x > 2) & (x < nx - 3) & (y > 2) & (y < ny - 3)
        table, x, y = table[on_chip], x[on_chip], y[on_chip]

        stars = {
            "x": np.asarray(x), "y": np.asarray(y),
            "ra": table["ra"].to_numpy(), "dec": table["dec"].to_numpy(),
            "source_id": table["source_id"].to_numpy(),
        }
        for band, column in CATALOG_COLUMNS.items():
            stars[band] = table[column].to_numpy()
        stars["B_V"] = stars["B"] - stars["V"]
        stars["source_id"] = stars["source_id"].astype("int64")

        faint = self._faint_tail(
            n_stars, bright_limit=float(np.nanmax(stars["V"])) if stars["V"].size else 18.0
        )
        for key in ("x", "y", "ra", "dec", "B", "V", "R", "I", "B_V", "source_id"):
            stars[key] = np.concatenate([stars[key], faint[key]])
        return stars

    def _faint_tail(self, n_stars, bright_limit=18.0, faint_limit=21.0, slope=0.35):
        """Synthetic stars below the catalogue limit, continuing the counts.

        Drawn from ``dN/dm ~ 10**(0.35 m)``, close enough to real Galactic counts
        that the detection limit falls where it would in practice. ``source_id``
        is -1, so a test can tell a catalogue star from an invented one.
        """
        ny, nx = self.shape
        area_scale = (nx * ny) / (1024.0 * 1024.0)
        count = max(int(n_stars * area_scale), 0)
        if count == 0:
            return {k: np.array([]) for k in
                    ("x", "y", "ra", "dec", "B", "V", "R", "I", "B_V", "source_id")}

        u = self.rng.random(count)
        lo, hi = 10 ** (slope * bright_limit), 10 ** (slope * faint_limit)
        v_mag = np.log10(lo + u * (hi - lo)) / slope

        x = self.rng.uniform(2, nx - 3, count)
        y = self.rng.uniform(2, ny - 3, count)
        b_v = self.rng.uniform(0.3, 1.5, count)
        v_r = 0.55 * b_v + 0.02 * self.rng.normal(size=count)
        ra, dec = self.wcs.all_pix2world(x, y, 0)
        return {
            "x": x, "y": y, "ra": np.asarray(ra), "dec": np.asarray(dec),
            "V": v_mag, "B": v_mag + b_v, "R": v_mag - v_r,
            "I": v_mag - v_r - 0.5 * b_v, "B_V": b_v,
            "source_id": np.full(count, -1, dtype="int64"),
        }

    def _build_stars_synthetic(self, n_stars):
        """A star field with a realistic magnitude distribution and colours.

        Counts rise toward faint magnitudes as ``10**(0.35 m)``, close enough to
        real Galactic counts that the detection limit falls where it would in
        practice. A handful of bright stars are forced in so that saturation and
        the zero point both have something to work with.
        """
        ny, nx = self.shape
        n_faint = max(n_stars - 12, 0)

        # Inverse-transform sampling of dN/dm ~ 10**(0.35 m) over [12, 20].
        u = self.rng.random(n_faint)
        m_lo, m_hi, slope = 12.0, 20.0, 0.35
        lo, hi = 10 ** (slope * m_lo), 10 ** (slope * m_hi)
        v_mag = np.log10(lo + u * (hi - lo)) / slope
        v_mag = np.concatenate([v_mag, self.rng.uniform(8.5, 12.0, min(12, n_stars))])

        x = self.rng.uniform(2, nx - 3, v_mag.size)
        y = self.rng.uniform(2, ny - 3, v_mag.size)

        # B-V spanning roughly A0 to K5, so a colour term is measurable.
        b_v = self.rng.uniform(-0.05, 1.35, v_mag.size)
        # V-R tracks B-V along the main sequence closely enough for this purpose.
        v_r = 0.55 * b_v + 0.02 * self.rng.normal(size=v_mag.size)

        ra, dec = self.wcs.all_pix2world(x, y, 0)
        return {
            "x": x, "y": y, "ra": ra, "dec": dec,
            "V": v_mag, "B": v_mag + b_v, "R": v_mag - v_r,
            "I": v_mag - v_r - 0.5 * b_v,
            "B_V": b_v,
        }

    def _build_galaxy(self):
        """One resolved Sersic galaxy near the centre.

        Its presence is what makes isophotal-vs-total magnitude, star/galaxy
        separation and (later) host-galaxy subtraction testable at all.
        """
        ny, nx = self.shape
        return {
            "x": nx * 0.5 + self.rng.uniform(-nx * 0.05, nx * 0.05),
            "y": ny * 0.5 + self.rng.uniform(-ny * 0.05, ny * 0.05),
            "V": 12.5,
            "B_V": 0.85,
            "r_eff": max(nx, ny) * 0.035,
            "n": 2.5,
            "ellip": 0.45,
            "theta": self.rng.uniform(0, np.pi),
        }

    # --- Rendering ------------------------------------------------------------
    def _magnitude(self, mags_or_value, band):
        """Per-band magnitude, from a V magnitude and a colour where needed."""
        if band in mags_or_value:
            return mags_or_value[band]
        return mags_or_value["V"]

    def flux_e_per_s(self, magnitude, band, airmass, transparency=1.0):
        """Above-atmosphere magnitude to detected e-/s, through the atmosphere."""
        extinguished = magnitude + EXTINCTION.get(band, 0.15) * airmass
        return transparency * 10 ** (-0.4 * (extinguished - ZERO_POINTS.get(band, 20.0)))

    def sky_e_per_s(self, band, airmass):
        """Sky background per pixel in e-/s, brightening with airmass."""
        per_arcsec2 = SKY_BRIGHTNESS.get(band, 19.0) - 0.2 * (airmass - 1.0)
        per_pixel = per_arcsec2 - 2.5 * np.log10(self.pixel_scale**2)
        return 10 ** (-0.4 * (per_pixel - ZERO_POINTS.get(band, 20.0)))

    def render(self, band, seeing_fwhm, airmass, transparency=1.0, dx=0.0, dy=0.0):
        """The scene in e-/s per pixel, at a given seeing, airmass and dither."""
        image = np.full(self.shape, self.sky_e_per_s(band, airmass), dtype=float)

        for index in range(self.stars["x"].size):
            x = self.stars["x"][index] + dx
            y = self.stars["y"][index] + dy
            flux = self.flux_e_per_s(self.stars[band][index], band, airmass, transparency)
            fwhm = self.local_fwhm(x, y, seeing_fwhm)
            _add_moffat(image, x, y, flux, fwhm)

        if self.galaxy is not None:
            self._add_galaxy(image, band, seeing_fwhm, airmass, transparency, dx, dy)
        return image

    def _add_galaxy(self, image, band, seeing_fwhm, airmass, transparency, dx, dy):
        """A Sersic profile, convolved with the PSF at the galaxy's position."""
        from astropy.convolution import Gaussian2DKernel, convolve_fft
        from astropy.modeling.models import Sersic2D

        galaxy = self.galaxy
        magnitude = galaxy["V"] if band == "V" else galaxy["V"] + {
            "B": galaxy["B_V"], "R": -0.5, "I": -0.9
        }.get(band, 0.0)
        total = self.flux_e_per_s(magnitude, band, airmass, transparency)

        ny, nx = self.shape
        y, x = np.mgrid[0:ny, 0:nx]
        model = Sersic2D(
            amplitude=1.0, r_eff=galaxy["r_eff"], n=galaxy["n"],
            x_0=galaxy["x"] + dx, y_0=galaxy["y"] + dy,
            ellip=galaxy["ellip"], theta=galaxy["theta"],
        )
        profile = np.asarray(model(x, y), dtype=float)
        profile[~np.isfinite(profile)] = 0.0
        if profile.sum() <= 0:
            return
        profile *= total / profile.sum()

        # Convolving the whole profile is affordable and avoids the artefacts of
        # rendering an unconvolved cusp.
        fwhm = self.local_fwhm(galaxy["x"], galaxy["y"], seeing_fwhm)
        kernel = Gaussian2DKernel(fwhm / 2.3548)
        image += convolve_fft(profile, kernel, normalize_kernel=True, allow_huge=True)

    def truth_table(self, bands):
        """The truth catalogue: one row per source, with total magnitudes."""
        rows = []
        for index in range(self.stars["x"].size):
            row = {
                "id": index + 1,
                "type": "star",
                "source_id": int(self.stars["source_id"][index])
                if "source_id" in self.stars else -1,
                "x": self.stars["x"][index],
                "y": self.stars["y"][index],
                "ra": self.stars["ra"][index],
                "dec": self.stars["dec"][index],
                "B_V": self.stars["B_V"][index],
            }
            for band in bands:
                row[f"mag_{band}"] = self.stars[band][index]
            rows.append(row)

        if self.galaxy is not None:
            row = {
                "id": len(rows) + 1, "type": "galaxy",
                "x": self.galaxy["x"], "y": self.galaxy["y"],
                "B_V": self.galaxy["B_V"],
            }
            ra, dec = self.wcs.all_pix2world(self.galaxy["x"], self.galaxy["y"], 0)
            row["ra"], row["dec"] = float(ra), float(dec)
            for band in bands:
                row[f"mag_{band}"] = self.galaxy["V"] + {
                    "B": self.galaxy["B_V"], "V": 0.0, "R": -0.5, "I": -0.9
                }.get(band, 0.0)
            rows.append(row)
        return rows


def _add_moffat(image, x, y, flux, fwhm, beta=MOFFAT_BETA, stamp_radius_fwhm=6.0):
    """Add one Moffat source to ``image`` in place, on a local stamp.

    The stamp is cut at 6xFWHM. A Moffat has genuinely extended wings -- that is
    the point of using one -- so this truncation loses a little flux; the tests
    account for it by comparing against :func:`enclosed_fraction` at the same
    radius rather than against the nominal total.
    """
    ny, nx = image.shape
    alpha = moffat_alpha(fwhm, beta)
    radius = max(int(np.ceil(stamp_radius_fwhm * fwhm)), 3)

    x0, x1 = max(int(x - radius), 0), min(int(x + radius) + 1, nx)
    y0, y1 = max(int(y - radius), 0), min(int(y + radius) + 1, ny)
    if x0 >= x1 or y0 >= y1:
        return

    yy, xx = np.mgrid[y0:y1, x0:x1]
    r2 = (xx - x) ** 2 + (yy - y) ** 2
    peak = flux * (beta - 1.0) / (np.pi * alpha**2)
    image[y0:y1, x0:x1] += peak * (1.0 + r2 / alpha**2) ** (-beta)
