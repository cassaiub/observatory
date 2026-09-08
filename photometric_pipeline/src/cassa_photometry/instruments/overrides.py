"""Layer user-configured detector constants over an instrument profile.

This is what makes ``detector:`` in a config file work. It builds a subclass of
whatever profile was chosen, so the result *is* an ``InstrumentProfile`` -- every
phase keeps working, and nothing downstream needs to know a user override is in
play.

The ordering is the part that matters scientifically. By default a configured
value fills in only where the frame's header is silent, preserving the
pipeline's rule that the frame is the first authority on its own data:

    header  ->  detector: block  ->  profile's own knowledge  ->  fallbacks

``override_header: true`` inverts the first two. That is occasionally the right
answer -- acquisition software does write wrong cards -- but it means a number
measured at the telescope is being discarded in favour of one typed into a YAML
file, so it is logged the first time it happens for each quantity.
"""

from cassa_photometry.logging_utils import get_logger

#: Quantities already warned about, so the log carries one line per run, not one
#: per frame. Keyed by (profile name, quantity).
_WARNED = set()


def _warn_once(profile_name, quantity, header_value, config_value):
    key = (profile_name, quantity)
    if key in _WARNED:
        return
    _WARNED.add(key)
    get_logger("cassa_photometry").warning(
        "detector.override_header is set: using the configured %s (%s) instead of "
        "the value in the frame header (%s).",
        quantity, config_value, header_value,
    )


def with_detector_overrides(profile, detector):
    """Return a profile that answers from ``detector`` as well as its own tables.

    Parameters
    ----------
    profile : InstrumentProfile
        The profile chosen by name.
    detector : DetectorOverride
        The parsed ``detector:`` config block.
    """
    base_cls = type(profile)
    override_header = bool(detector.override_header)

    def _resolve(quantity, header, configured, header_value, profile_value):
        """Apply the order of authority to one quantity.

            header  ->  detector: block  ->  profile's own knowledge

        The middle position is the whole point. Asking the profile's ``get_*``
        for the header value would not work: a profile like ``Cassa8InchProfile``
        answers from its own detector curve when the header is silent, so it
        never returns ``None`` and a configured value would be shadowed by a
        table the user was trying to override. Hence the ``*_from_header`` seam.
        """
        if header_value is not None and not (override_header and configured is not None):
            return header_value
        if configured is not None:
            if header_value is not None and header_value != configured:
                _warn_once(base_cls.__name__, quantity, header_value, configured)
            return configured
        return profile_value

    class ConfiguredProfile(base_cls):
        """``{base}`` plus the detector constants from the config file."""

        #: The config block this profile was built from, for logging and tests.
        DETECTOR_OVERRIDE = detector

        @property
        def name(self):
            return f"{super().name} + configured detector"

        def get_gain(self, header):
            return _resolve(
                "gain", header, detector.gain,
                self.gain_from_header(header), super().get_gain(header),
            )

        def get_read_noise(self, header):
            return _resolve(
                "read noise", header, detector.read_noise,
                self.read_noise_from_header(header), super().get_read_noise(header),
            )

        def get_saturation(self, header):
            return _resolve(
                "saturation", header, detector.saturation_adu,
                self.saturation_from_header(header), super().get_saturation(header),
            )

        def pixel_scale_arcsec(self):
            # The profile's optics are the fallback; a configured scale wins
            # here because a user who types one knows their own telescope.
            if detector.pixel_scale_arcsec is not None:
                return detector.pixel_scale_arcsec
            return super().pixel_scale_arcsec()

        def flat_proxies(self):
            proxies = dict(super().flat_proxies())
            proxies.update(
                {str(k).upper(): str(v).upper() for k, v in detector.flat_proxies.items()}
            )
            return proxies

    # Filter tables are merged rather than replaced, so a user adding one odd
    # filter name does not lose the standard ones.
    if detector.filter_map:
        ConfiguredProfile.FILTER_MAP = {
            **base_cls.FILTER_MAP,
            **{str(k).upper(): str(v) for k, v in detector.filter_map.items()},
        }
    if detector.science_bands:
        ConfiguredProfile.SCIENCE_BANDS = {
            **base_cls.SCIENCE_BANDS,
            **{str(k).upper(): str(v).upper() for k, v in detector.science_bands.items()},
        }
    if detector.uncalibrated_filters:
        ConfiguredProfile.UNCALIBRATED_FILTERS = set(base_cls.UNCALIBRATED_FILTERS) | {
            str(f).upper() for f in detector.uncalibrated_filters
        }

    ConfiguredProfile.__doc__ = ConfiguredProfile.__doc__.replace("{base}", base_cls.__name__)
    ConfiguredProfile.__name__ = f"Configured{base_cls.__name__}"
    return ConfiguredProfile()
