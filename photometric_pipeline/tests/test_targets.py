"""Resolving what is being observed.

The frame's own header is the authority, because the observer knew what they
were pointing at and a card cannot be lost the way a side file can. Everything
here is about precedence and about failing softly: an unknown target is a field,
never an error.
"""

import numpy as np
import pytest
from conftest import header

from cassa_photometry.targets import (
    FIELD,
    GALAXY,
    SUPERNOVA,
    VARIABLE,
    Target,
    TargetRegistry,
    normalise_type,
)


def _registry(tmp_path=None, text=None):
    if text is None:
        return TargetRegistry()
    path = tmp_path / "targets.yaml"
    path.write_text(text)
    return TargetRegistry(str(path))


# --- The header is the authority ----------------------------------------------

def test_a_frame_declaring_its_target_needs_no_targets_file():
    target = _registry().resolve(header(
        TARGNAME="SN2026abc", OBJTYPE="supernova",
        TARGRA="22 37 06.21", TARGDEC="+34 24 51.08",
        HOSTGAL="NGC7331", DISCDATE="2026-09-01",
    ))
    assert target.type == SUPERNOVA
    assert target.name == "SN2026abc"
    assert target.host == "NGC7331"
    assert target.discovery == "2026-09-01"
    assert target.ra == pytest.approx(339.2759, abs=1e-3)
    assert target.dec == pytest.approx(34.4142, abs=1e-3)
    assert target.source == "header"


def test_object_alone_is_a_named_field():
    target = _registry().resolve(header(OBJECT="NGC7331"))
    assert target.type == FIELD
    assert target.name == "NGC7331"


def test_a_frame_with_nothing_is_a_field_not_an_error():
    """Today's behaviour, and the safe one: a missing target must never stop a
    reduction."""
    target = _registry().resolve(header())
    assert target.type == FIELD
    assert target.is_field


def test_an_objtype_without_a_name_falls_back_to_a_field():
    target = _registry().resolve(header(OBJTYPE="supernova", OBJECT="NGC7331"))
    assert target.type == FIELD


def test_an_unrecognised_objtype_is_a_field_and_says_so(pipeline_logs):
    target = _registry().resolve(header(TARGNAME="X", OBJTYPE="quasar-ish"))
    assert target.type == FIELD


@pytest.mark.parametrize("raw, expected", [
    ("SN", SUPERNOVA), ("supernova", SUPERNOVA), ("Transient", SUPERNOVA),
    ("var", VARIABLE), ("Variable Star", VARIABLE),
    ("galaxy", GALAXY), ("extended", GALAXY),
    ("field", FIELD), ("STAR", "star"),
])
def test_objtype_spellings_are_normalised(raw, expected):
    assert normalise_type(raw) == expected


def test_an_unknown_spelling_is_none_rather_than_a_guess():
    assert normalise_type("nonsense") is None
    assert normalise_type(None) is None


# --- The targets file overrides and supplements --------------------------------

def test_a_targets_file_overrides_the_header(tmp_path):
    """Its real job: correcting a refined position after discovery."""
    registry = _registry(tmp_path, """
targets:
  - {name: SN2026abc, type: supernova, ra: 339.5, dec: 34.5}
""")
    target = registry.resolve(header(TARGNAME="SN2026abc", OBJTYPE="supernova",
                                     TARGRA=339.0, TARGDEC=34.0))
    assert target.ra == pytest.approx(339.5)
    assert "targets.yaml" in target.source


def test_a_targets_file_supplies_what_the_telescope_cannot_know(tmp_path):
    registry = _registry(tmp_path, """
targets:
  - {name: V0532Cyg, type: variable, comparison: [GSC 3162-1155, GSC 3162-0902]}
""")
    target = registry.resolve(header(OBJECT="V0532Cyg"))
    assert target.type == VARIABLE
    assert len(target.comparison) == 2


def test_a_targets_file_can_retrofit_frames_that_predate_the_cards(tmp_path):
    registry = _registry(tmp_path, """
targets:
  - {name: NGC7331, type: galaxy}
""")
    target = registry.resolve(header(OBJECT="NGC 7331"))
    assert target.type == GALAXY


def test_names_match_regardless_of_spacing_and_case(tmp_path):
    registry = _registry(tmp_path, "targets:\n  - {name: ngc_7331, type: galaxy}\n")
    assert registry.resolve(header(OBJECT="NGC 7331")).type == GALAXY


def test_a_target_not_in_the_file_still_resolves_from_its_header(tmp_path):
    registry = _registry(tmp_path, "targets:\n  - {name: SomethingElse, type: galaxy}\n")
    target = registry.resolve(header(TARGNAME="SN2026abc", OBJTYPE="supernova"))
    assert target.type == SUPERNOVA


def test_a_missing_targets_file_warns_and_carries_on(pipeline_logs):
    registry = TargetRegistry("/nowhere/targets.yaml")
    assert registry.resolve(header(OBJECT="NGC7331")).type == FIELD
    assert "not found" in pipeline_logs.text


def test_a_malformed_targets_file_does_not_stop_a_reduction(tmp_path, pipeline_logs):
    path = tmp_path / "targets.yaml"
    path.write_text("targets: [ this is not: valid: yaml\n")
    registry = TargetRegistry(str(path))
    assert registry.resolve(header(OBJECT="NGC7331")).type == FIELD


# --- Description --------------------------------------------------------------

def test_a_target_describes_itself_for_the_log():
    assert "field" in Target(name="NGC7331").describe()
    assert "supernova SN2026abc at" in Target(
        name="SN2026abc", type=SUPERNOVA, ra=339.0, dec=34.0, source="header"
    ).describe()


# --- The empirical PSF --------------------------------------------------------

def test_an_epsf_is_built_from_a_frames_own_stars():
    from cassa_photometry.psf import build_epsf, estimate_fwhm
    from cassa_photometry.simulate.scene import _add_moffat

    rng = np.random.default_rng(0)
    image = rng.normal(100.0, 5.0, (400, 400))
    for x, y in rng.uniform(40, 360, (25, 2)):
        _add_moffat(image, x, y, 60000.0, 4.2)

    model, info = build_epsf(image, estimate_fwhm(image, fwhm_guess=4.0, threshold=5.0))
    assert model is not None
    assert info["n_stars"] >= 4


def test_a_star_poor_field_degrades_rather_than_crashing():
    """A sparse high-latitude field is a normal thing to observe; every caller
    falls back to an analytic Gaussian of the measured FWHM."""
    from cassa_photometry.psf import build_epsf

    model, info = build_epsf(np.zeros((64, 64)), {"positions": [], "fwhm_px": 4.0})
    assert model is None
    assert "need at least 4" in info["reason"]

    model, info = build_epsf(np.zeros((64, 64)), None)
    assert model is None and info["reason"]


def test_the_psf_module_is_importable_from_its_old_location():
    """Phase 4 and existing notebooks import it from there."""
    from cassa_photometry.phase4_diagnostics.psf import estimate_fwhm as old
    from cassa_photometry.psf import estimate_fwhm as new

    assert old is new
