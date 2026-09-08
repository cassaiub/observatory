"""Instrument-profile registry: the one place a profile is chosen by name.

Every phase resolves its profile through :func:`get_profile` from
``config.instrument``, so a run is reduced with the profile the user asked for
rather than one hardcoded at import time.

There are three ways to describe an observing setup, and a new one should reach
for the first that suffices:

1. **A ``detector:`` block in the config file.** No code at all. Covers the
   common case of a camera that does not write its own gain and read noise.
   See :class:`~cassa_photometry.config.DetectorOverride`.
2. **A profile class in a local file**, named by ``instrument_module:`` in the
   config as ``path/to/file.py:ClassName``. Needed when the answer depends on
   the header -- a CMOS conversion-gain curve, a multi-amplifier layout.
3. **A profile shipped by another installed package**, which registers itself
   under the ``cassa_photometry.instruments`` entry-point group.

Only setups the CASSA Observatory itself operates are built in, so the registry
stays a statement about this observatory rather than a directory of everyone's
telescopes.
"""

from cassa_photometry.instruments.base import InstrumentProfile
from cassa_photometry.instruments.cassa import Cassa8InchProfile

#: Config name -> profile class, for setups shipped with the pipeline.
PROFILES = {
    # Standard FITS keywords only; correct for any setup that writes them, and
    # the right starting point for a new rig when paired with a ``detector:``
    # config block.
    "generic": InstrumentProfile,
    "cassa8": Cassa8InchProfile,
}

#: Entry-point group third-party packages register profiles under.
ENTRY_POINT_GROUP = "cassa_photometry.instruments"

DEFAULT_PROFILE = "generic"


def _plugin_profiles():
    """Profile classes contributed by other installed packages.

    Discovered lazily and defensively: a broken third-party plugin must not stop
    the pipeline from reducing data with a built-in profile.
    """
    from importlib.metadata import entry_points

    found = {}
    try:
        points = entry_points(group=ENTRY_POINT_GROUP)
    except Exception:  # pragma: no cover - importlib metadata problems
        return found
    for point in points:
        try:
            found[point.name.strip().lower()] = point.load()
        except Exception as exc:  # pragma: no cover - depends on third-party code
            import logging

            logging.getLogger(__name__).warning(
                "Instrument plugin %r failed to load: %s", point.name, exc
            )
    return found


def _load_from_path(spec):
    """Load ``path/to/file.py:ClassName`` and return the class.

    Lets a user keep their own profile beside their data, with no packaging and
    no edit to this repository -- which is what keeps a local clone rebasable.
    """
    import importlib.util
    import os

    path, _, class_name = spec.partition(":")
    path = os.path.abspath(os.path.expanduser(path))
    if not class_name:
        raise KeyError(
            f"instrument_module {spec!r} must be written as 'path/to/file.py:ClassName'."
        )
    if not os.path.isfile(path):
        raise KeyError(f"instrument_module file not found: {path}")

    module_spec = importlib.util.spec_from_file_location("cassa_user_profile", path)
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    try:
        return getattr(module, class_name)
    except AttributeError:
        raise KeyError(f"{path} defines no class named {class_name!r}.") from None


def available_profiles():
    """Sorted registry names, for CLI help text and error messages."""
    return sorted(set(PROFILES) | set(_plugin_profiles()))


def get_profile(name=None, config=None):
    """Return an instrument profile instance.

    Parameters
    ----------
    name : str, optional
        Registry name, or a ``path.py:ClassName`` spec. ``None`` yields the
        ``generic`` profile.
    config : PipelineConfig, optional
        When given, ``config.instrument_module`` can name a local profile file
        and ``config.detector`` is layered on top of whichever profile is
        chosen.
    """
    spec = None
    if config is not None and getattr(config, "instrument_module", None):
        spec = config.instrument_module
    elif name and ".py:" in name:
        spec = name

    if spec:
        profile = _load_from_path(spec)()
    else:
        key = (name or DEFAULT_PROFILE).strip().lower()
        classes = {**PROFILES, **_plugin_profiles()}
        try:
            profile = classes[key]()
        except KeyError:
            raise KeyError(
                f"Unknown instrument profile {name!r}. "
                f"Available: {', '.join(sorted(classes))}. "
                f"To describe a setup with no profile, use a 'detector:' config "
                f"block or 'instrument_module: my_profile.py:MyProfile' "
                f"(see docs/CUSTOMIZING.md)."
            ) from None

    detector = getattr(config, "detector", None) if config is not None else None
    if detector is not None and not detector.is_empty():
        from cassa_photometry.instruments.overrides import with_detector_overrides

        profile = with_detector_overrides(profile, detector)
    return profile
