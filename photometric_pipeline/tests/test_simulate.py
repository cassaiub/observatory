"""The simulator, and the properties that make it worth testing against.

A simulator is only useful if it is faithful in the ways the pipeline is being
checked on, and reproducible so a failure can be re-run. These assert both --
including the two properties that a naive simulator gets wrong and that would
silently invalidate every downstream photometric test.
"""

import numpy as np
import pytest
from astropy.io import fits

from cassa_photometry.simulate import PRESETS, run
from cassa_photometry.simulate.detector import DetectorModel
from cassa_photometry.simulate.presets import CASSA_8INCH, delivered_fwhm_arcsec
from cassa_photometry.simulate.scene import MOFFAT_BETA, Scene, enclosed_fraction


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """One tiny dataset, shared: generating it is the slow part."""
    outdir = tmp_path_factory.mktemp("sim")
    return run(preset="tiny", outdir=str(outdir), seed=4242)


# --- Reproducibility ----------------------------------------------------------

def test_the_same_seed_gives_the_same_pixels(tmp_path):
    """The dataset is distributed as a seed, so this is the distribution."""
    import glob
    import os

    first = run(preset="tiny", outdir=str(tmp_path / "a"), seed=7)
    second = run(preset="tiny", outdir=str(tmp_path / "b"), seed=7)
    assert first["n_frames"] == second["n_frames"]

    for path in sorted(glob.glob(os.path.join(first["raw_dir"], "*.fits"))):
        twin = os.path.join(second["raw_dir"], os.path.basename(path))
        assert np.array_equal(fits.getdata(path), fits.getdata(twin)), path


def test_a_different_seed_gives_different_pixels(tmp_path):
    a = run(preset="tiny", outdir=str(tmp_path / "a"), seed=1)
    b = run(preset="tiny", outdir=str(tmp_path / "b"), seed=2)
    import glob
    import os

    name = os.path.basename(sorted(glob.glob(os.path.join(a["raw_dir"], "NGC*.fits")))[0])
    assert not np.array_equal(
        fits.getdata(os.path.join(a["raw_dir"], name)),
        fits.getdata(os.path.join(b["raw_dir"], name)),
    )


# --- Faithfulness -------------------------------------------------------------

def test_the_psf_is_a_moffat_with_real_wings():
    """A Gaussian PSF has ~zero aperture correction, which is why a Gaussian
    test suite passes while the aperture correction is missing."""
    fwhm = 4.2
    moffat = enclosed_fraction(2.0 * fwhm, fwhm, beta=MOFFAT_BETA)
    assert 0.92 < moffat < 0.95
    # ...and the correction it implies is the ~0.074 mag the pipeline must find.
    assert -2.5 * np.log10(moffat) == pytest.approx(0.074, abs=0.005)


def test_the_psf_broadens_off_axis():
    """A constant PSF cannot show that a scalar aperture correction is wrong."""
    scene = Scene((512, 512), 0.598, 339.266875, 34.415778,
                  rng=np.random.default_rng(0), field="ngc7331", galaxy=False)
    centre = scene.local_fwhm(255, 255, 4.0)
    corner = scene.local_fwhm(0, 0, 4.0)
    assert centre == pytest.approx(4.0, abs=0.01)
    assert corner > 1.3 * centre


def test_the_star_field_is_real_so_a_solver_can_use_it():
    """Random star positions would leave phase 2 and the zero point untested."""
    scene = Scene((1024, 1024), 0.598, 339.266875, 34.415778,
                  rng=np.random.default_rng(0), field="ngc7331", n_stars=200)
    real = scene.stars["source_id"] > 0
    assert real.sum() > 50, "too few catalogue stars for a plate solve"
    # The faint tail continues the counts below Gaia's synthetic-photometry
    # limit; without it the frame is sparser than the real sky.
    assert (~real).sum() > 0
    assert scene.stars["V"][real].min() < 14.0


def test_rendered_flux_matches_the_requested_magnitude():
    """The truth catalogue is only truth if the pixels agree with it."""
    scene = Scene((256, 256), 0.598, 339.266875, 34.415778,
                  rng=np.random.default_rng(0), field=None, n_stars=1, galaxy=False)
    scene.stars = {k: np.array(v) for k, v in {
        "x": [128.0], "y": [128.0], "ra": [0.0], "dec": [0.0],
        "V": [15.0], "B": [15.5], "R": [14.7], "I": [14.4], "B_V": [0.5],
    }.items()}

    expected = scene.flux_e_per_s(15.0, "V", airmass=1.0)
    image = scene.render("V", seeing_fwhm=4.0, airmass=1.0)
    measured = float(np.sum(image - scene.sky_e_per_s("V", 1.0)))
    # The stamp is cut at 6xFWHM, so a little of the Moffat's wings is lost.
    assert measured == pytest.approx(expected, rel=0.01)


def test_extinction_dims_sources_with_airmass():
    scene = Scene((64, 64), 0.598, 0.0, 0.0, rng=np.random.default_rng(0),
                  field=None, n_stars=1, galaxy=False)
    high = scene.flux_e_per_s(15.0, "V", airmass=1.0)
    low = scene.flux_e_per_s(15.0, "V", airmass=2.0)
    assert low < high
    assert -2.5 * np.log10(low / high) == pytest.approx(0.15, abs=1e-6)


def test_seeing_is_not_the_delivered_psf():
    """Sizing apertures from the DIMM number under-sizes them by ~30%."""
    delivered = delivered_fwhm_arcsec(2.0, 1.5, "V", guide_rms=0.42, optics_fwhm=1.0)
    assert delivered == pytest.approx(2.87, abs=0.05)
    assert delivered > 1.4 * 2.0


# --- The detector model -------------------------------------------------------

def test_non_linearity_grows_with_signal():
    """A brightness-dependent error tilts a magnitude scale; a constant one
    would merely offset it, and would not matter."""
    detector = DetectorModel((16, 16), rng=np.random.default_rng(0), nonlinearity=0.02)
    low = detector._apply_nonlinearity(np.array([5000.0]))[0]
    high = detector._apply_nonlinearity(np.array([50000.0]))[0]
    assert (5000.0 - low) / 5000.0 < (50000.0 - high) / 50000.0


def test_frames_never_exceed_full_well():
    detector = DetectorModel((32, 32), rng=np.random.default_rng(0))
    saturated = detector.expose_flat(1e6, 10.0)
    assert saturated.max() <= detector.full_well_adu


# --- What is written out ------------------------------------------------------

def test_frames_carry_the_cards_the_pipeline_reads(dataset):
    import glob
    import os

    science = sorted(glob.glob(os.path.join(dataset["raw_dir"], "NGC*.fits")))[0]
    header = fits.getheader(science)

    for card in ("EGAIN", "READNOIS", "SATURATE", "SECPIX", "XPIXSZ", "FOCALLEN",
                 "FILTER", "EXPTIME", "IMAGETYP", "OBJECT", "OBJCTRA", "OBJCTDEC",
                 "DATE-OBS", "MJD-OBS", "AIRMASS", "SITELAT", "SITELONG",
                 "SEEING", "SEEINGWL", "SEEINGER", "GUIDERMS",
                 "TARGNAME", "OBJTYPE", "TARGRA", "TARGDEC"):
        assert card in header, f"missing {card}"

    # Raw means raw: phase 1 refuses frames that claim prior calibration.
    assert header["CALSTAT"] == ""
    assert header["BUNIT"] == "ADU"
    # And no WCS -- phase 2 has to solve the field, as it does with real data.
    assert "CD1_1" not in header and "CRVAL1" not in header
    # A simulated frame must never be mistakable for real data.
    assert header["SIMULATD"] is True


def test_calibration_frames_are_typed_and_flats_say_what_kind(dataset):
    import glob
    import os

    kinds = {}
    for path in glob.glob(os.path.join(dataset["raw_dir"], "*.fits")):
        header = fits.getheader(path)
        kinds.setdefault(header["IMAGETYP"], []).append(header)

    assert {"BIAS", "DARK", "FLAT", "LIGHT"} <= set(kinds)
    assert all(h["EXPTIME"] == 0.0 for h in kinds["BIAS"])
    # FLATTYPE distinguishes a sky flat (corrects illumination) from a dome or
    # panel flat (does not), which changes what the flat is entitled to remove.
    assert all(h["FLATTYPE"] == "sky" for h in kinds["FLAT"])


def test_the_truth_files_describe_every_source_and_frame(dataset):
    import pandas as pd

    sources = pd.read_csv(dataset["truth_sources"])
    frames = pd.read_csv(dataset["truth_frames"])

    assert len(sources) == dataset["n_sources"]
    assert {"id", "type", "x", "y", "ra", "dec", "mag_V"} <= set(sources.columns)
    assert len(frames) == dataset["n_science"]
    assert {"mjd", "band", "airmass", "delivered_fwhm_px", "transparency",
            "true_zp"} <= set(frames.columns)
    assert frames["delivered_fwhm_px"].between(1.0, 30.0).all()


def test_twilight_flats_fade_through_the_sequence(dataset):
    """A median across raw flats would track the fading sky, not the pixel
    response -- which is why phase 1 normalises each flat before combining."""
    import glob
    import os

    levels = [
        float(np.median(fits.getdata(p)))
        for p in sorted(glob.glob(os.path.join(dataset["raw_dir"], "flat_*.fits")))
    ]
    assert levels[0] > levels[-1] * 1.1


def test_an_unknown_preset_names_the_alternatives(tmp_path):
    with pytest.raises(KeyError) as exc:
        run(preset="nope", outdir=str(tmp_path))
    for name in PRESETS:
        assert name in str(exc.value)


def test_the_workshop_preset_is_large_enough_to_plate_solve():
    """A preset that only solved with an unusual index set would be a trap:
    the shipped 4203-4206 indexes need a field of roughly 10 arcmin or more."""
    ny, nx = PRESETS["workshop"]["shape"]
    scale = CASSA_8INCH["pixel_scale"]
    assert min(nx, ny) * scale / 60.0 > 10.0
