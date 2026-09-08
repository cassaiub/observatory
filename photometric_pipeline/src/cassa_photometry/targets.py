"""What is being observed, and therefore how it should be measured.

A star, a galaxy and a supernova need different photometry, and the pipeline
cannot choose sensibly without knowing which it is looking at. That knowledge
exists at the telescope -- the observer pointed at something on purpose -- so it
belongs in the frame's header, where it is one card and cannot be lost.
Reconstructing it afterwards from a side file is possible but lossy.

Resolution order, most authoritative first:

1. **The frame's own header** -- ``TARGNAME``, ``OBJTYPE``, ``TARGRA``/``TARGDEC``,
   ``HOSTGAL``, ``DISCDATE``.
2. **A ``targets.yaml``**, which overrides and supplements it. Its real jobs are
   retrofitting archival frames that predate those cards, correcting a refined
   position after discovery, and holding what the telescope cannot know -- an
   explicit comparison-star list, for instance.
3. **``OBJECT`` alone**, treated as a named field.
4. **Nothing**, which is a field.

A missing or unrecognised target is never an error. It falls through to field
mode, which is exactly today's behaviour.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

from cassa_photometry.logging_utils import get_logger

#: Target types the pipeline routes on.
STAR = "star"
VARIABLE = "variable"
SUPERNOVA = "supernova"
GALAXY = "galaxy"
FIELD = "field"

KNOWN_TYPES = (STAR, VARIABLE, SUPERNOVA, GALAXY, FIELD)

#: Header cards that describe the target.
TARGET_CARDS = ("TARGNAME", "OBJTYPE", "TARGRA", "TARGDEC", "HOSTGAL", "DISCDATE")


@dataclass
class Target:
    """One observed target."""

    name: str = ""
    type: str = FIELD
    ra: float | None = None
    dec: float | None = None
    #: Host galaxy, for supernova host modelling and template selection.
    host: str | None = None
    #: Discovery date, which excludes post-discovery epochs from a template
    #: stack -- a deep stack containing the supernova would self-subtract it.
    discovery: str | None = None
    #: Explicit comparison stars for differential photometry. The telescope
    #: cannot know these, so they only ever come from a targets file.
    comparison: list = field(default_factory=list)
    #: Where this target came from, for the log and the product header.
    source: str = "default"

    @property
    def is_field(self):
        return self.type == FIELD

    @property
    def has_position(self):
        return self.ra is not None and self.dec is not None

    def describe(self):
        if self.is_field:
            return f"field {self.name or '(unnamed)'}"
        position = (f" at {self.ra:.5f} {self.dec:+.5f}" if self.has_position else "")
        return f"{self.type} {self.name}{position} [{self.source}]"


def normalise_type(value):
    """A header's ``OBJTYPE`` mapped onto a known type, or None."""
    if value in (None, ""):
        return None
    text = str(value).strip().lower()
    aliases = {
        "sn": SUPERNOVA, "supernova": SUPERNOVA, "transient": SUPERNOVA,
        "var": VARIABLE, "variable": VARIABLE, "variable star": VARIABLE,
        "star": STAR, "point": STAR, "standard": STAR,
        "galaxy": GALAXY, "gal": GALAXY, "extended": GALAXY,
        "field": FIELD, "survey": FIELD,
    }
    return aliases.get(text)


class TargetRegistry:
    """Resolves the target of a frame."""

    def __init__(self, targets_file=None, logger=None):
        self.logger = logger or get_logger("cassa_photometry")
        self._overrides = _load_targets_file(targets_file, self.logger)

    def resolve(self, header):
        """The :class:`Target` a frame describes."""
        target = self._from_header(header)
        override = self._match_override(target, header)
        if override is not None:
            target = _merge(target, override)
        if target.type not in KNOWN_TYPES:
            self.logger.warning(
                "Unknown OBJTYPE %r; treating this as a field. Known types: %s.",
                target.type, ", ".join(KNOWN_TYPES),
            )
            target.type = FIELD
        return target

    # -- resolution steps ------------------------------------------------------
    def _from_header(self, header):
        if header is None:
            return Target(source="default")

        name = _text(header.get("TARGNAME"))
        declared = normalise_type(header.get("OBJTYPE"))
        if declared and name:
            return Target(
                name=name, type=declared,
                ra=_degrees(header.get("TARGRA"), is_ra=True),
                dec=_degrees(header.get("TARGDEC"), is_ra=False),
                host=_text(header.get("HOSTGAL")) or None,
                discovery=_text(header.get("DISCDATE")) or None,
                source="header",
            )
        if declared and not name:
            self.logger.info(
                "OBJTYPE is %r but TARGNAME is missing; treating this as a field. "
                "Both cards are needed to route a target.", header.get("OBJTYPE"),
            )
        # A named field is still worth recording, even with no type.
        object_name = _text(header.get("OBJECT"))
        return Target(name=object_name, type=FIELD,
                      source="OBJECT" if object_name else "default")

    def _match_override(self, target, header):
        """The targets-file entry for this frame, by target name then OBJECT."""
        if not self._overrides:
            return None
        for candidate in (target.name, _text(header.get("OBJECT")) if header else ""):
            key = _key(candidate)
            if key and key in self._overrides:
                return self._overrides[key]
        return None


def _merge(base, override):
    """A targets-file entry layered over what the header said."""
    merged = Target(
        name=override.name or base.name,
        type=override.type if override.type != FIELD else base.type,
        ra=override.ra if override.ra is not None else base.ra,
        dec=override.dec if override.dec is not None else base.dec,
        host=override.host or base.host,
        discovery=override.discovery or base.discovery,
        comparison=list(override.comparison or base.comparison),
        source="targets.yaml" if base.source == "default" else f"{base.source}+targets.yaml",
    )
    return merged


def _load_targets_file(path, logger):
    """Parse a ``targets.yaml``, keyed by normalised name."""
    if not path:
        return {}
    if not os.path.exists(path):
        logger.warning("Targets file not found: %s. Falling back to the headers.", path)
        return {}
    try:
        import yaml

        with open(path) as handle:
            data = yaml.safe_load(handle) or {}
    except Exception as exc:
        logger.warning("Could not read %s (%s); falling back to the headers.", path, exc)
        return {}

    entries = data.get("targets", data if isinstance(data, list) else [])
    resolved = {}
    for entry in entries or []:
        if not isinstance(entry, dict) or not entry.get("name"):
            continue
        target = Target(
            name=str(entry["name"]),
            type=normalise_type(entry.get("type")) or FIELD,
            ra=_degrees(entry.get("ra"), is_ra=True),
            dec=_degrees(entry.get("dec"), is_ra=False),
            host=entry.get("host"),
            discovery=str(entry["discovery"]) if entry.get("discovery") else None,
            comparison=list(entry.get("comparison") or []),
            source="targets.yaml",
        )
        resolved[_key(target.name)] = target
    logger.info("Loaded %d target(s) from %s.", len(resolved), path)
    return resolved


def _key(name):
    """Names are matched case- and space-insensitively, as OBJECT cards vary."""
    return str(name or "").strip().upper().replace(" ", "").replace("_", "")


def _text(value):
    return str(value).strip() if value not in (None, "") else ""


def _degrees(value, is_ra):
    """Degrees from a decimal or sexagesimal coordinate, or None."""
    if value in (None, ""):
        return None
    from cassa_photometry.phase2_integration.wcs import _as_degrees

    return _as_degrees(value, is_ra=is_ra)
