"""The declared steps, the registry and the code must agree.

Three things have to line up for a toggle to mean anything:

1. the ``PhaseNSteps`` dataclass declares it, so it can be set from YAML;
2. ``cassa_photometry.steps`` declares it, so it has a place in the order and a
   dependency contract;
3. the phase actually consults its plan, so setting it changes what runs.

Any one of those missing gives a setting that validates and does nothing --
which is the defect this whole design exists to remove, and which sat unnoticed
in seven toggles because ``test_config.py`` only ever proved that the YAML
parses. These tests fail if the three drift apart.
"""

import re
from pathlib import Path

import pytest

from cassa_photometry import config as config_module
from cassa_photometry import steps as steps_module
from cassa_photometry.config import (
    Phase1Steps,
    Phase2Steps,
    Phase3Steps,
    Phase4Steps,
    StepToggles,
)

SRC = Path(config_module.__file__).parent

#: (config attribute, toggles class, implementing package).
STEP_CLASSES = [
    ("phase1", Phase1Steps, "phase1_calibration"),
    ("phase2", Phase2Steps, "phase2_integration"),
    ("phase3", Phase3Steps, "phase3_photometry"),
    ("phase4", Phase4Steps, "phase4_diagnostics"),
]


def _all_step_names():
    return [
        (phase, name, package)
        for phase, cls, package in STEP_CLASSES
        for name in cls.step_names()
    ]


def _phase_source(package):
    """Every line of a phase's own implementation.

    ``config.py`` is excluded because it *declares* the toggles; matching there
    would make the search trivially self-satisfying.
    """
    chunks = []
    for path in sorted((SRC / package).rglob("*.py")):
        if path.name == "config.py" or "__pycache__" in path.parts:
            continue
        chunks.append(path.read_text())
    return "\n".join(chunks)


def _consumes(phase, step, source):
    """Whether ``source`` gates anything on ``step``.

    Two shapes count, because the codebase legitimately uses both:

    1. A literal gate -- ``steps.enabled("cosmic_rays")``.
    2. A plan-driven phase, which asks ``cassa_photometry.steps`` for the
       resolved order and dispatches on the names it gets back. The step name
       never appears inside a call, so a literal-only search would report a
       correctly wired toggle as dead.

    For (2) the requirement is still two-sided: the phase must resolve *its own*
    plan and mention the step name somewhere as a string. A name that appears
    nowhere in the phase is still caught.
    """
    literal = re.compile(rf"""steps\.enabled\(\s*["']{re.escape(step)}["']\s*\)""")
    if literal.search(source):
        return True
    plan_driven = re.compile(
        rf"""resolve(?:d_names|_plan)\(\s*["']{re.escape(phase)}["']""")
    named = re.compile(rf"""["']{re.escape(step)}["']""")
    return bool(plan_driven.search(source) and named.search(source))


# --- The three-way agreement --------------------------------------------------

@pytest.mark.parametrize(("phase", "step", "package"), _all_step_names(),
                         ids=lambda v: str(v))
def test_every_declared_step_is_consumed_by_its_own_phase(phase, step, package):
    """A declared toggle with no consumer is a lie told to the user."""
    assert _consumes(phase, step, _phase_source(package)), (
        f"{phase}.steps.{step} is declared in config.py but nothing in "
        f"{package}/ consumes it. Setting it to false validates, logs nothing "
        f"and runs the step anyway. Either wire it up or remove the "
        f"declaration -- see docs/step-plan.md."
    )


@pytest.mark.parametrize(("phase", "step", "package"), _all_step_names(),
                         ids=lambda v: str(v))
def test_every_declared_step_is_in_the_registry(phase, step, package):
    """Without a registry entry a step has no place in the order and no
    dependency contract, so it could never be reordered or validated."""
    registered = {s.name for s in steps_module.registry_for(phase)}
    assert step in registered, (
        f"{phase}.steps.{step} is declared in config.py but missing from "
        f"cassa_photometry.steps.{phase.upper()}."
    )


@pytest.mark.parametrize(("phase", "cls", "package"), STEP_CLASSES,
                         ids=lambda v: str(v))
def test_every_registry_step_is_configurable(phase, cls, package):
    """The reverse: a registry entry nobody can switch off is not a step, it is
    just code with a name."""
    declared = set(cls.step_names())
    for step in steps_module.registry_for(phase):
        assert step.name in declared, (
            f"cassa_photometry.steps.{phase.upper()} declares {step.name!r} "
            f"but {cls.__name__} has no such field, so it cannot be configured."
        )


# --- The guard's own integrity ------------------------------------------------

def test_the_guard_rejects_a_toggle_nothing_gates_on():
    """A check that cannot fail is worse than no check: it reads as coverage."""
    assert not _consumes("phase1", "never_wired", 'steps.enabled("something")')
    # A name mentioned but never gated on does not count as consumed.
    assert not _consumes("phase1", "never_wired", 'x = "never_wired"')
    # A plan-driven phase must resolve its OWN plan, not another's.
    assert not _consumes(
        "phase1", "never_wired",
        'resolved_names("phase2", t)\nx = "never_wired"')
    assert _consumes(
        "phase1", "never_wired",
        'resolved_names("phase1", t)\nx = "never_wired"')


def test_the_toggle_classes_are_all_steptoggles():
    for _, cls, _ in STEP_CLASSES:
        assert issubclass(cls, StepToggles)


def test_control_fields_are_not_mistaken_for_steps():
    """`exclude`, `order` and `custom` configure the plan; they are not steps,
    and must not be reported as skipped ones."""
    for _, cls, _ in STEP_CLASSES:
        for control in StepToggles._CONTROL:
            assert control not in cls.step_names()


def test_no_phase_declares_zero_steps():
    for phase, cls, _ in STEP_CLASSES:
        assert cls.step_names(), f"{phase} declares no steps at all"


def test_each_phase_package_exists():
    """A renamed package must fail loudly rather than making the guard vacuous."""
    for _, _, package in STEP_CLASSES:
        assert (SRC / package).is_dir(), f"{package}/ not found under {SRC}"
