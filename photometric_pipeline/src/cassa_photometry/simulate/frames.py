"""Turn a scene and a detector into a night's worth of FITS frames.

The headers are as much the deliverable as the pixels. A simulated dataset whose
frames omit the cards the pipeline reads would validate the arithmetic while
hiding every metadata bug, so these frames carry the full acquisition spec --
including the time-domain and seeing cards that nothing reads *yet* but that the
epoch and PSF work depends on.

One card is deliberately absent: raw frames carry **no WCS**. The pointing is
given as ``OBJCTRA``/``OBJCTDEC`` only, so phase 2 has to solve the field the way
it does with real data. That works because the star field is real.
"""

import os

import numpy as np
from astropy.io import fits
from astropy.time import Time

#: Written into every simulated frame, so a frame can never be mistaken for real
#: data -- by a person or by a reduction.
SIMULATED_CARD = ("SIMULATD", True, "Synthetic frame from cassa-simulate")


def _base_header(observation, detector, instrument):
    """The cards every frame carries, whatever its type."""
    header = fits.Header()
    header["SIMULATD"] = (True, "Synthetic frame from cassa-simulate")
    header["SIMVERS"] = (observation["sim_version"], "cassa-simulate version")
    header["SIMSEED"] = (observation["seed"], "RNG seed; the dataset is reproducible")

    header["TELESCOP"] = (instrument["telescope"], "Telescope")
    header["INSTRUME"] = (instrument["camera"], "Camera")
    header["FOCALLEN"] = (instrument["focal_length_mm"], "[mm] Focal length")
    header["XPIXSZ"] = (instrument["pixel_size_um"], "[um] Pixel size, after binning")
    header["YPIXSZ"] = (instrument["pixel_size_um"], "[um] Pixel size, after binning")
    header["SECPIX"] = (instrument["pixel_scale"], "[arcsec/pixel] Plate scale")
    header["XBINNING"] = (1, "Binning factor, X")
    header["YBINNING"] = (1, "Binning factor, Y")

    # Detector constants the frame states about itself. A real camera that writes
    # these makes an instrument profile unnecessary, which is the point of
    # including them.
    header["EGAIN"] = (detector.gain, "[e-/ADU] System gain")
    header["READNOIS"] = (detector.read_noise_e, "[e-] Read noise")
    header["GAIN"] = (instrument["gain_setting"], "Unitless gain setting (NOT e-/ADU)")
    header["READOUTM"] = (instrument["read_mode"], "Readout mode")
    header["SATURATE"] = (detector.full_well_adu, "[ADU] Saturation level")
    header["CCD-TEMP"] = (instrument["ccd_temp"], "[C] Sensor temperature")
    header["SET-TEMP"] = (instrument["ccd_temp"], "[C] Cooler set point")
    header["BUNIT"] = ("ADU", "Raw counts")
    header["CALSTAT"] = ("", "No calibration applied")

    header["SITELAT"] = (instrument["site_lat"], "[deg] Site latitude")
    header["SITELONG"] = (instrument["site_long"], "[deg] Site longitude, East positive")
    header["TIMESYS"] = ("UTC", "Time scale")
    return header


def _add_time(header, when, exposure):
    """``DATE-OBS``/``MJD-OBS``/``DATE-END`` for one exposure."""
    from astropy.time import TimeDelta

    start = Time(when, scale="utc")
    header["DATE-OBS"] = (start.isot, "UTC at exposure start")
    header["MJD-OBS"] = (float(start.mjd), "MJD at exposure start")
    header["DATE-END"] = (
        (start + TimeDelta(exposure, format="sec")).isot, "UTC at exposure end"
    )
    header["EXPTIME"] = (float(exposure), "[s] Exposure time")
    header["EXPOSURE"] = (float(exposure), "[s] Exposure time")
    return header


def write_science(path, data, scene, observation, detector, instrument, target):
    """A science frame, with pointing, target, and observing-condition cards."""
    header = _base_header(observation, detector, instrument)
    _add_time(header, observation["when"], observation["exposure"])

    header["IMAGETYP"] = ("LIGHT", "Frame type")
    header["FILTER"] = (observation["band"], "Filter")
    header["OBJECT"] = (target["object"], "Field name")

    # Where the telescope pointed, in degrees. No WCS: phase 2 solves the field.
    header["OBJCTRA"] = (observation["pointing_ra"], "[deg] Commanded RA")
    header["OBJCTDEC"] = (observation["pointing_dec"], "[deg] Commanded Dec")
    header["AIRMASS"] = (observation["airmass"], "Airmass at mid-exposure")

    # Target cards: what is being observed, known at the telescope.
    header["TARGNAME"] = (target["name"], "Canonical target designation")
    header["OBJTYPE"] = (target["type"], "star|variable|supernova|galaxy|field")
    header["TARGRA"] = (target["ra"], "[deg] Catalogued target RA")
    header["TARGDEC"] = (target["dec"], "[deg] Catalogued target Dec")

    # Seeing, from the DIMM. Zenith-corrected at SEEINGWL -- NOT the delivered
    # image FWHM, which is larger by the airmass, wavelength, guiding and optics
    # terms. Recording the convention is what keeps the two from being confused.
    header["SEEING"] = (observation["dimm_seeing"], "[arcsec] DIMM zenith seeing")
    header["SEEINGWL"] = (500.0, "[nm] Wavelength SEEING refers to")
    header["SEEINGER"] = (observation["dimm_seeing_err"], "[arcsec] DIMM seeing error")
    header["GUIDERMS"] = (observation["guide_rms"], "[arcsec] RMS guiding error")

    # Truth, for tests. Real frames do not carry these.
    header["SIMFWHM"] = (observation["fwhm_px"], "[pixel] TRUTH: delivered FWHM at centre")
    header["SIMTRANS"] = (observation["transparency"], "TRUTH: relative transparency")
    header["SIMZP"] = (observation["true_zp"], "TRUTH: instrumental zero point")

    _write(path, data, header)
    return path


def write_calibration(path, data, kind, observation, detector, instrument, band=None,
                      flat_type=None):
    """A bias, dark or flat frame."""
    header = _base_header(observation, detector, instrument)
    _add_time(header, observation["when"], observation["exposure"])
    header["IMAGETYP"] = (kind.upper(), "Frame type")
    header["OBJECT"] = (kind.upper(), "Calibration frame")
    if band is not None:
        header["FILTER"] = (band, "Filter")
    if flat_type is not None:
        # Sky flats correct illumination as well as pixel response; dome and
        # panel flats do not. Recording which is what lets a reduction say
        # whether a large-scale gradient has been corrected or absorbed.
        header["FLATTYPE"] = (flat_type, "sky|dome|panel")
    _write(path, data, header)
    return path


def _write(path, data, header):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fits.PrimaryHDU(np.asarray(data, dtype=np.uint16), header).writeto(path, overwrite=True)


def airmass_from_altitude(altitude_deg):
    """Plane-parallel airmass with the usual Hardie-style refraction correction."""
    z = np.radians(90.0 - np.asarray(altitude_deg, dtype=float))
    sec_z = 1.0 / np.cos(z)
    return sec_z - 0.0018167 * (sec_z - 1) - 0.002875 * (sec_z - 1) ** 2 \
        - 0.0008083 * (sec_z - 1) ** 3
