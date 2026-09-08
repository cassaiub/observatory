"""Photometric systems: which system each band is on, and what that costs.

Two distinct errors used to live here, both silent and both larger than the
pipeline's quoted zero-point uncertainty.

**Vega versus AB.** The flux-calibrated image converts magnitudes to Janskys
with the AB relation, ``F = 3631 Jy x 10**(-m/2.5)``, for every band. But B and
V are calibrated against APASS *Johnson* magnitudes, which are Vega-based. A
Vega magnitude fed through the AB relation is wrong by the band's AB offset:
**B by about 8%, V by about 2%**.

**Johnson versus Sloan.** `FILTER='R'` is a Johnson-Cousins R filter, and it was
being calibrated against APASS Sloan r'. Measured on real Gaia DR3 synthetic
photometry for stars in the NGC 7331 field, **Johnson R minus SDSS r' is
-0.21 +/- 0.036 mag**: a fifth of a magnitude of systematic zero-point error,
plus colour-dependent scatter, in two of the five supported bands.

The fix for both is to say plainly which system a band is on and to calibrate
against a reference in that same system. Gaia DR3 synthetic photometry supplies
Johnson-Cousins *and* SDSS magnitudes for the same stars, which is why it is the
preferred anchor.
"""

#: AB minus Vega magnitude, per band. A Vega magnitude plus this offset is the
#: AB magnitude, which is what the Jy conversion actually requires.
#: Values from Blanton & Roweis (2007) / Bessell (1990) synthesis; they are
#: properties of the filter curves, not of this instrument.
AB_MINUS_VEGA = {
    "U": 0.79,
    "B": -0.09,
    "V": 0.02,
    "R": 0.21,
    "I": 0.45,
    # Sloan bands are defined on the AB system, so the offset is zero.
    "G": 0.0,
}

#: Which system each science band's *filter* is on.
BAND_SYSTEM = {
    "U": "Vega", "B": "Vega", "V": "Vega", "R": "Vega", "I": "Vega",
    "G": "AB",
}

#: Effective wavelength per band in nm, recorded so a Jy value states the
#: wavelength its monochromatic flux density refers to.
BAND_WAVELENGTH_NM = {
    "U": 365.0, "B": 445.0, "V": 551.0, "R": 658.0, "I": 806.0, "G": 464.0,
}


def system_of(band):
    """``"Vega"`` or ``"AB"`` for a science band; ``"AB"`` when unknown."""
    return BAND_SYSTEM.get(str(band).upper(), "AB")


def ab_offset(band):
    """Magnitudes to add to convert this band's calibrated magnitude to AB.

    Zero for a band already on AB. Applying it is what makes the Jy conversion
    correct rather than off by up to 8%.
    """
    return float(AB_MINUS_VEGA.get(str(band).upper(), 0.0))


def effective_wavelength_nm(band):
    """Effective wavelength of a band in nm, or None."""
    return BAND_WAVELENGTH_NM.get(str(band).upper())


def describe(band):
    """A one-line summary for a header comment or a log line."""
    band = str(band).upper()
    system = system_of(band)
    offset = ab_offset(band)
    if offset:
        return f"{band}: {system} system, AB-Vega = {offset:+.2f} mag"
    return f"{band}: {system} system"
