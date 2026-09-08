"""Simulated CASSA 8-inch data with known truth.

Generates a complete observing run -- bias, dark, flat and science frames --
from a real star field, through a forward model that is the deliberate inverse
of what the pipeline undoes. Every source's true magnitude, every frame's true
seeing, transparency and zero point are written out alongside, so a reduction can
be scored against truth rather than against plausibility.

    cassa-simulate --preset workshop --out sim/

Why generate rather than distribute: the dataset is a seed, not a download. A
student runs one command and gets the same frames as everyone else, and the
tests can afford to regenerate a small version on every run.

Two properties are what make this worth testing against, and both are choices
that a naive simulator gets wrong:

* **The star field is real** (Gaia DR3 positions with synthetic Johnson-Cousins
  magnitudes), so ``solve-field`` genuinely solves the frames and the
  reference-catalog cross-match genuinely finds the stars. A random star field
  would leave phase 2 and the zero point untested.
* **The PSF is a Moffat that broadens off-axis.** A Gaussian PSF has an aperture
  correction of essentially zero, so a Gaussian test suite passes while the
  aperture correction is missing -- which is exactly what happened here.
"""

import os

import numpy as np

from cassa_photometry.logging_utils import get_logger
from cassa_photometry.simulate.detector import DetectorModel
from cassa_photometry.simulate.frames import write_calibration, write_science
from cassa_photometry.simulate.presets import (
    CASSA_8INCH,
    PRESETS,
    TARGET,
    NightModel,
    delivered_fwhm_arcsec,
)
from cassa_photometry.simulate.scene import ZERO_POINTS, Scene

__all__ = ["run", "PRESETS", "Scene", "DetectorModel"]

SIM_VERSION = "1.0"

#: Night 1 starts here; later nights are spaced a day apart. The frames
#: deliberately straddle UT midnight, because binning epochs on the UT date
#: string splits a single observing night into two.
FIRST_NIGHT_START = "2026-09-02T21:30:00"


def run(preset="workshop", outdir="sim", seed=20260902, config=None, logger=None,
        overwrite=True):
    """Generate a simulated dataset. Returns a summary dict.

    Parameters
    ----------
    preset : str
        One of :data:`PRESETS` -- ``tiny``, ``workshop`` or ``multinight``.
    outdir : str
        Written to ``<outdir>/raw`` plus three truth files at the top level.
    seed : int
        Everything derives from this, so a dataset is reproducible.
    """
    from astropy.time import Time, TimeDelta

    logger = logger or get_logger("cassa_simulate")
    if preset not in PRESETS:
        raise KeyError(f"Unknown preset {preset!r}. Available: {', '.join(PRESETS)}.")
    spec = PRESETS[preset]

    rng = np.random.default_rng(seed)
    raw_dir = os.path.join(outdir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    instrument = dict(CASSA_8INCH)
    detector = DetectorModel(spec["shape"], rng=rng)
    scene = Scene(
        spec["shape"], instrument["pixel_scale"], TARGET["ra"], TARGET["dec"],
        rng=rng, field=spec["field"], galaxy=spec["galaxy"],
        n_stars=spec.get("n_faint", 400),
    )

    logger.info("Simulating preset %r -> %s", preset, os.path.abspath(outdir))
    logger.info(
        "  %s, %d x %d at %.3f\"/px, %d stars on chip, seed %d",
        instrument["telescope"], spec["shape"][1], spec["shape"][0],
        instrument["pixel_scale"], scene.stars["x"].size, seed,
    )

    frame_rows = []
    written = []

    for night in range(spec["nights"]):
        night_start = Time(FIRST_NIGHT_START) + TimeDelta(night, format="jd")
        model = NightModel(rng, night_index=night)
        stamp = night_start.isot[:10].replace("-", "")

        written += _write_calibrations(
            raw_dir, spec, detector, instrument, night_start, stamp, night, rng, seed
        )

        # Science frames, cycling through the filters the way an observer would.
        n_total = spec["n_science"] * len(spec["bands"])
        index = 0
        for repeat in range(spec["n_science"]):
            for band in spec["bands"]:
                conditions = model.step((index + 0.5) / max(n_total, 1))
                when = night_start + TimeDelta(
                    (index + 1) * (spec["exposure"] + 20.0), format="sec"
                )
                row = _write_science_frame(
                    raw_dir, scene, detector, instrument, spec, band, conditions,
                    when, stamp, repeat, night, rng, seed, logger,
                )
                frame_rows.append(row)
                written.append(row["file"])
                index += 1

    truth = _write_truth(outdir, scene, spec, frame_rows, instrument, detector, seed, preset)
    logger.info("Wrote %d frames and %d truth sources.", len(written), truth["n_sources"])
    logger.info("Next:  cassa-run -i %s -o %s --instrument cassa8",
                raw_dir, os.path.join(outdir, "work"))
    return {
        "outdir": os.path.abspath(outdir),
        "raw_dir": os.path.abspath(raw_dir),
        "n_frames": len(written),
        "n_science": len(frame_rows),
        **truth,
    }


def _write_calibrations(raw_dir, spec, detector, instrument, night_start, stamp,
                        night, rng, seed):
    """Bias, dark and flat sequences for one night."""
    from astropy.time import TimeDelta

    written = []
    observation = {"sim_version": SIM_VERSION, "seed": seed,
                   "when": night_start.isot, "exposure": 0.0}

    for i in range(spec["n_bias"]):
        observation["when"] = (night_start + TimeDelta(i * 5.0, format="sec")).isot
        observation["exposure"] = 0.0
        path = os.path.join(raw_dir, f"bias_{stamp}_n{night}_{i + 1:03d}.fits")
        written.append(write_calibration(
            path, detector.expose_bias(), "BIAS", observation, detector, instrument
        ))

    for i in range(spec["n_dark"]):
        observation["when"] = (night_start + TimeDelta(60.0 + i * 70.0, format="sec")).isot
        observation["exposure"] = spec["exposure"]
        path = os.path.join(raw_dir, f"dark_{stamp}_n{night}_{i + 1:03d}.fits")
        written.append(write_calibration(
            path, detector.expose_dark(spec["exposure"]), "DARK",
            observation, detector, instrument
        ))

    # Twilight flats fade as they are taken. That is what makes normalising each
    # flat by its own median before combining the right order of operations, so
    # the simulator reproduces the fade rather than assuming it away.
    for band in spec["bands"]:
        level = 30000.0
        for i in range(spec["n_flat"]):
            level *= rng.uniform(0.88, 0.94)
            observation["when"] = (
                night_start - TimeDelta(1800.0 - i * 30.0, format="sec")
            ).isot
            observation["exposure"] = 3.0
            path = os.path.join(raw_dir, f"flat_{band}_{stamp}_n{night}_{i + 1:03d}.fits")
            written.append(write_calibration(
                path, detector.expose_flat(level, 3.0), "FLAT",
                observation, detector, instrument, band=band, flat_type="sky"
            ))
    return written


def _write_science_frame(raw_dir, scene, detector, instrument, spec, band, conditions,
                         when, stamp, repeat, night, rng, seed, logger):
    """One science exposure, and the truth row describing it."""
    fwhm_arcsec = delivered_fwhm_arcsec(
        conditions["dimm_seeing"], conditions["airmass"], band,
        conditions["guide_rms"], instrument["optics_fwhm_arcsec"],
    )
    fwhm_px = fwhm_arcsec / instrument["pixel_scale"]

    # Dithering, so registration has something to do and hot pixels do not land
    # on the same sources in every frame.
    dx, dy = rng.normal(0.0, 3.0, 2)

    photon_rate = scene.render(
        band, fwhm_px, conditions["airmass"],
        transparency=conditions["transparency"], dx=dx, dy=dy,
    )
    data = detector.expose_science(photon_rate, spec["exposure"])

    # The zero point this frame is actually on: the instrumental zero point,
    # dimmed by extinction and transparency, for a 1-second exposure.
    from cassa_photometry.simulate.scene import EXTINCTION

    true_zp = (
        ZERO_POINTS.get(band, 20.0)
        - EXTINCTION.get(band, 0.15) * conditions["airmass"]
        + 2.5 * np.log10(max(conditions["transparency"], 1e-6))
    )

    observation = {
        "sim_version": SIM_VERSION, "seed": seed, "when": when.isot,
        "exposure": spec["exposure"], "band": band,
        "pointing_ra": TARGET["ra"] + dx * instrument["pixel_scale"] / 3600.0,
        "pointing_dec": TARGET["dec"] + dy * instrument["pixel_scale"] / 3600.0,
        "airmass": conditions["airmass"],
        "dimm_seeing": conditions["dimm_seeing"],
        "dimm_seeing_err": conditions["dimm_seeing_err"],
        "guide_rms": conditions["guide_rms"],
        "transparency": conditions["transparency"],
        "fwhm_px": fwhm_px,
        "true_zp": true_zp,
    }
    name = f"{TARGET['object']}_{band}_{stamp}_n{night}_{repeat + 1:03d}.fits"
    path = os.path.join(raw_dir, name)
    write_science(path, data, scene, observation, detector, instrument, TARGET)

    return {
        "file": path, "night": night, "band": band, "date_obs": when.isot,
        "mjd": float(when.mjd), "exposure": spec["exposure"],
        "airmass": conditions["airmass"], "dimm_seeing": conditions["dimm_seeing"],
        "delivered_fwhm_arcsec": fwhm_arcsec, "delivered_fwhm_px": fwhm_px,
        "transparency": conditions["transparency"], "true_zp": true_zp,
        "dither_dx": float(dx), "dither_dy": float(dy),
    }


def _write_truth(outdir, scene, spec, frame_rows, instrument, detector, seed, preset):
    """The three truth files a reduction is scored against."""
    import pandas as pd
    import yaml

    sources = pd.DataFrame(scene.truth_table(spec["bands"]))
    sources_path = os.path.join(outdir, "truth_sources.csv")
    sources.to_csv(sources_path, index=False)

    frames_path = os.path.join(outdir, "truth_frames.csv")
    pd.DataFrame(frame_rows).to_csv(frames_path, index=False)

    config_path = os.path.join(outdir, "truth_config.yaml")
    with open(config_path, "w") as handle:
        yaml.safe_dump(
            {
                "preset": preset,
                "seed": seed,
                "sim_version": SIM_VERSION,
                "instrument": {k: v for k, v in instrument.items()},
                "detector": {
                    "gain_e_per_adu": detector.gain,
                    "read_noise_e": detector.read_noise_e,
                    "bias_level_adu": detector.bias_level,
                    "dark_rate_e_per_s": detector.dark_rate,
                    "full_well_adu": detector.full_well_adu,
                    "nonlinearity_at_full_well": detector.nonlinearity,
                },
                "photometry": {
                    "zero_points": ZERO_POINTS,
                    "moffat_beta": 2.5,
                    "coma_strength": scene.coma_strength,
                },
                "note": (
                    "Truth for a simulated dataset. mag_* columns in "
                    "truth_sources.csv are TOTAL magnitudes above the atmosphere, "
                    "in Johnson-Cousins."
                ),
            },
            handle,
            sort_keys=False,
        )
    return {
        "n_sources": len(sources),
        "truth_sources": sources_path,
        "truth_frames": frames_path,
        "truth_config": config_path,
    }
