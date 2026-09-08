"""The step registry: what each phase does, in what order, and what may move.

A phase is a named, ordered list of steps. Each step declares what it
``requires`` and what it ``provides``, which is what makes a *custom* order
checkable rather than merely accepted.

The rule, in one sentence: **if step S requires token T and some included step P
provides T, then P must run before S.** Two consequences follow, and both are
deliberate:

* **Excluding a step is always legal.** If nothing provides T, nothing requires
  it either -- reducing frames that were flat-fielded elsewhere is a real
  workflow, not an error to be caught.
* **Misordering included steps is never legal.** Flat-fielding before bias
  subtraction produces a plausible, silently wrong image, which is exactly the
  class of failure this pipeline already spends `CALSKIP`, `CRVETO` and the
  pre-calibrated-frame guard on preventing.

Tokens encode *physics*, not taste. Where two orders are both defensible --
cosmic rays before or after flat fielding, background subtraction before or
after registration, the aperture correction before or after the zero point --
no token links them, and the plan resolver permits either. Over-constraining
here would be as wrong as under-constraining: it would refuse reductions that
are perfectly sound.
"""

from __future__ import annotations

from dataclasses import dataclass, field


class StepPlanError(ValueError):
    """A step plan asked for something the pipeline cannot honour."""


@dataclass(frozen=True)
class Step:
    """One named operation in a phase."""

    name: str
    #: Tokens that must already exist. Only constrains steps that are included.
    requires: tuple = ()
    #: Tokens this step makes true.
    provides: tuple = ()
    #: Set for steps whose position is fixed by the phase's control flow rather
    #: than by physics -- they cannot be moved even though nothing forbids it on
    #: dependency grounds. Saying so is better than silently ignoring an order.
    movable: bool = True
    #: Loaded callable, for a user-supplied step. None for built-ins.
    function: object = field(default=None, compare=False)


# --- Phase 1: instrument signature removal ------------------------------------
#
# The order below is the canonical one. `linearity` is applied in raw ADU before
# anything else touches the pixel values, `overscan` changes the array shape so
# nothing may precede it that assumes the final geometry, and bias -> dark ->
# flat is the standard ISR chain: each correction assumes the previous one has
# happened. `cosmic_rays` deliberately requires only a trimmed frame, so it may
# be placed before or after flat fielding -- both are defensible and the choice
# is the user's.
PHASE1 = (
    Step("linearity", provides=("linearised",)),
    Step("overscan", requires=("linearised",), provides=("trimmed",)),
    Step("bad_pixel_mask", provides=("bpm",)),
    Step("bias", requires=("trimmed",), provides=("bias_subtracted",)),
    Step("dark", requires=("bias_subtracted",), provides=("dark_subtracted",)),
    Step("flat", requires=("dark_subtracted",), provides=("flat_fielded",)),
    Step("cosmic_rays", requires=("trimmed", "bpm"), provides=("cr_cleaned",)),
    Step("measure_fwhm", requires=("flat_fielded",), provides=("frame_fwhm",)),
)

# --- Phase 2: integration -----------------------------------------------------
#
# `align` and `subtract_background` are independent: a background model can be
# fitted before or after registration, and neither is obviously right. Both must
# precede `stack`, which consumes their output.
PHASE2 = (
    Step("subtract_background", provides=("background_removed",)),
    Step("align", provides=("registered",)),
    Step("stack", requires=("registered", "background_removed"),
         provides=("stacked",)),
    Step("solve_wcs", requires=("stacked",), provides=("wcs",)),
    # Measured on the master and written into its header, so it necessarily
    # happens between the combine and the write. There is no seam to move it to.
    Step("measure_fwhm", requires=("stacked",), provides=("master_fwhm",),
         movable=False),
    Step("visual_qa", requires=("stacked", "wcs"), provides=("qa_pdf",)),
)

# --- Phase 3: photometry ------------------------------------------------------
#
# `aperture_correction` is intentionally NOT required by `zero_point`: measuring
# the correction from the curve of growth and applying it to a combined zero
# point are both reasonable orders, and the pipeline supports either.
# `flux_calibration` genuinely cannot precede the zero point -- there is nothing
# to scale the image by -- and `classification` is defined as MAG_PSF minus
# MAG_AUTO, so it cannot precede the PSF photometry that produces MAG_PSF.
#
# Three of these run *inside* another step rather than beside it: the aperture
# correction is measured during the zero point, and PSF photometry and
# classification during catalog generation. They are switchable but not
# movable -- there is no seam to move them to, and pretending otherwise would
# accept an `order` the phase then ignores.
PHASE3 = (
    Step("aperture_correction", provides=("apcor",), movable=False),
    Step("zero_point", provides=("zero_point",)),
    Step("flux_calibration", requires=("zero_point",), provides=("fluxcal",)),
    Step("catalog", provides=("catalog",)),
    Step("psf_photometry", requires=("catalog",), provides=("psf_mag",),
         movable=False),
    Step("classification", requires=("catalog", "psf_mag"),
         provides=("classes",), movable=False),
)

# --- Phase 4: diagnostics -----------------------------------------------------
#
# The four stages inspect four different products and share no state, so every
# order is valid. They are listed in pipeline order because that is the order a
# reader of the report expects, not because anything requires it.
PHASE4 = (
    Step("stage_0_raw", provides=("stage_0",)),
    Step("stage_1_calibrated", provides=("stage_1",)),
    Step("stage_2_master", provides=("stage_2",)),
    Step("stage_3_photometry", provides=("stage_3",)),
)

REGISTRY = {
    "phase1": PHASE1,
    "phase2": PHASE2,
    "phase3": PHASE3,
    "phase4": PHASE4,
}


def registry_for(phase):
    """The declared steps of a phase, in canonical order."""
    try:
        return REGISTRY[phase]
    except KeyError:
        raise StepPlanError(
            f"Unknown phase {phase!r}. Known: {', '.join(sorted(REGISTRY))}."
        ) from None


def load_custom_step(name, spec):
    """Load a user-supplied step from ``path/to/file.py:callable``.

    Deliberately the same convention as ``instrument_module``, so there is one
    thing to learn. The callable may carry ``requires`` / ``provides``
    attributes; without them the step is unconstrained and sits wherever the
    ``order`` puts it.
    """
    import importlib.util
    import os

    if not isinstance(spec, str) or ":" not in spec:
        raise StepPlanError(
            f"steps.custom.{name}: expected 'path/to/file.py:callable', "
            f"got {spec!r}."
        )
    path, _, attribute = spec.rpartition(":")
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise StepPlanError(f"steps.custom.{name}: no such file: {path}")

    module_spec = importlib.util.spec_from_file_location(
        f"cassa_custom_step_{name}", path)
    if module_spec is None or module_spec.loader is None:
        raise StepPlanError(f"steps.custom.{name}: cannot import {path}")
    module = importlib.util.module_from_spec(module_spec)
    try:
        module_spec.loader.exec_module(module)
    except Exception as exc:
        raise StepPlanError(f"steps.custom.{name}: {path} failed to import: {exc}") from exc

    function = getattr(module, attribute, None)
    if function is None:
        raise StepPlanError(f"steps.custom.{name}: {path} defines no {attribute!r}.")
    if not callable(function):
        raise StepPlanError(f"steps.custom.{name}: {attribute!r} is not callable.")

    return Step(
        name=name,
        requires=tuple(getattr(function, "requires", ()) or ()),
        provides=tuple(getattr(function, "provides", ()) or ()) + (name,),
        function=function,
    )


def _validate_names(phase, label, names, known):
    unknown = [n for n in names if n not in known]
    if unknown:
        raise StepPlanError(
            f"{phase}.steps.{label}: unknown step(s) "
            f"{', '.join(repr(n) for n in unknown)}. "
            f"Known here: {', '.join(sorted(known))}."
        )


def _check_fixed_positions(phase, declared, included):
    """Refuse to move a step whose position is set by the phase's control flow.

    Some steps are nested inside another -- phase 3's star/galaxy classification
    happens *during* catalog generation, not as a sibling of it. Such a step can
    be switched off, but not resequenced, because there is no seam to move it
    to. Refusing is the honest response: silently accepting an ``order`` the
    phase then ignores would be exactly the "config that validates and does
    nothing" defect this whole design exists to remove.
    """
    canonical = {step.name: i for i, step in enumerate(declared)}
    actual = {step.name: i for i, step in enumerate(included)}

    for step in included:
        if step.movable:
            continue
        for other in included:
            if other.name == step.name:
                continue
            # A custom step has no canonical position, so it cannot "move" a
            # built-in relative to itself. Comparing against one would make an
            # immovable step forbid every insertion point but the very end,
            # which would leave custom steps all but unusable in phase 3.
            if other.function is not None or other.name not in canonical:
                continue
            canonical_before = canonical[other.name] < canonical[step.name]
            actual_before = actual[other.name] < actual[step.name]
            if canonical_before != actual_before:
                raise StepPlanError(
                    f"{phase}.steps.order moves {step.name!r} relative to "
                    f"{other.name!r}, but {step.name!r} runs inside another "
                    f"step rather than as a sibling of it, so its position is "
                    f"fixed by this phase's control flow. It can be switched "
                    f"off, but not resequenced."
                )


def _check_dependencies(phase, ordered):
    """Refuse an order that puts a step before something it depends on."""
    provided_by = {}
    for step in ordered:
        for token in step.provides:
            provided_by.setdefault(token, step.name)

    seen = set()
    for step in ordered:
        for token in step.requires:
            # Nothing included provides this token, so nothing constrains us:
            # the capability is simply absent, which excluding a step is
            # *supposed* to be able to do.
            if token not in provided_by:
                continue
            if token not in seen:
                raise StepPlanError(
                    f"{phase}.steps.order: {step.name!r} is placed before "
                    f"{provided_by[token]!r}, but it needs what that step "
                    f"produces ({token!r}). Move {provided_by[token]!r} "
                    f"earlier, or exclude it if you do not want it at all."
                )
        seen.update(step.provides)


def resolve_plan(phase, toggles):
    """Resolve a phase's configuration into an ordered list of :class:`Step`.

    Applies, in order: the boolean toggles and ``exclude`` (which step runs),
    ``custom`` (what else exists), and ``order`` (in what sequence), then
    validates the result against the dependency graph.
    """
    declared = list(registry_for(phase))
    known = {step.name for step in declared}

    custom = dict(getattr(toggles, "custom", None) or {})
    for name, spec in custom.items():
        if name in known:
            raise StepPlanError(
                f"{phase}.steps.custom.{name}: a built-in step is already "
                f"called {name!r}; choose another name."
            )
        declared.append(load_custom_step(name, spec))
        known.add(name)

    exclude = list(getattr(toggles, "exclude", None) or [])
    _validate_names(phase, "exclude", exclude, known)

    included = [
        step for step in declared
        if step.name not in exclude and toggles.enabled(step.name)
    ]

    order = list(getattr(toggles, "order", None) or [])
    if order:
        _validate_names(phase, "order", order, known)
        wanted = [name for name in order if name not in exclude]
        included_names = [step.name for step in included]

        missing = [n for n in included_names if n not in wanted]
        extra = [n for n in wanted if n not in included_names]
        if missing or extra:
            problems = []
            if missing:
                problems.append(
                    f"does not mention {', '.join(repr(n) for n in missing)}")
            if extra:
                problems.append(
                    f"mentions {', '.join(repr(n) for n in extra)}, which "
                    f"is switched off")
            raise StepPlanError(
                f"{phase}.steps.order must list exactly the steps that run, "
                f"once each -- it {'; and '.join(problems)}. Run a command with "
                f"--show-plan to print the list to start from."
            )
        if len(set(wanted)) != len(wanted):
            raise StepPlanError(
                f"{phase}.steps.order repeats a step; each may appear once."
            )

        by_name = {step.name: step for step in included}
        included = [by_name[name] for name in wanted]

        _check_fixed_positions(phase, declared, included)

    _check_dependencies(phase, included)
    return included


def resolved_names(phase, toggles):
    """Just the names, in resolved order. The common case."""
    return [step.name for step in resolve_plan(phase, toggles)]


def describe_plan(phase, toggles):
    """A human-readable plan, for ``--show-plan``."""
    declared = registry_for(phase)
    try:
        running = resolved_names(phase, toggles)
    except StepPlanError as exc:
        return f"{phase}: INVALID -- {exc}"

    lines = [f"{phase}:"]
    for position, name in enumerate(running, 1):
        builtin = any(s.name == name for s in declared)
        lines.append(f"  {position:>2}. {name}" + ("" if builtin else "   [custom]"))
    for step in declared:
        if step.name not in running:
            lines.append(f"   -  {step.name}   [off]")
    return "\n".join(lines)
