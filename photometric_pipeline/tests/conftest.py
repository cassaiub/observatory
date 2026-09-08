"""Shared fixtures and helpers for the test suite.

The important one is :class:`FixtureProfile`. Many tests need *a* profile that
reports known detector constants -- they are testing calibration arithmetic, not
any particular telescope. That role used to be played by the iTelescope profile,
which coupled a dozen phase-1 tests to a setup the observatory does not operate.
A purpose-built fixture says what it is.
"""

import numpy as np
import pytest
from astropy.io import fits

from cassa_photometry.instruments.base import InstrumentProfile

#: Detector constants the fixture profile reports when a header is silent. These
#: are the values the phase-1 tests were written against, so keeping them means
#: the calibration assertions did not have to change.
FIXTURE_GAIN = 1.0
FIXTURE_READ_NOISE = 10.0


class FixtureProfile(InstrumentProfile):
    """A profile with known constants, standing in for a real camera.

    Behaves exactly like ``generic`` except that it answers instead of returning
    ``None`` when the header says nothing -- which is what a profile for a real
    camera does, and what the calibration tests need in order to have a defined
    error budget.
    """

    GAIN_E_PER_ADU = FIXTURE_GAIN
    READ_NOISE_E = FIXTURE_READ_NOISE

    @property
    def name(self):
        return "Test Fixture Instrument"

    def get_gain(self, header):
        gain = super().get_gain(header)
        return self.GAIN_E_PER_ADU if gain is None else gain

    def get_read_noise(self, header):
        read_noise = super().get_read_noise(header)
        return self.READ_NOISE_E if read_noise is None else read_noise


@pytest.fixture
def profile():
    """A profile with known detector constants."""
    return FixtureProfile()


def header(**cards):
    """A FITS header from keyword arguments."""
    return fits.Header(cards)


def write_raw(path, data, **cards):
    """Write a single-HDU raw frame with the given header cards."""
    hdu = fits.PrimaryHDU(np.asarray(data, dtype=np.float32), header(**cards))
    hdu.writeto(str(path), overwrite=True)
    return str(path)


@pytest.fixture
def pipeline_logs(caplog):
    """Make ``caplog`` see the pipeline's own loggers.

    ``get_logger`` sets ``propagate = False`` so phase output is not duplicated
    through the root logger, and it re-applies that on every call. caplog
    attaches at the root, so it would otherwise capture nothing and a test
    asserting on a warning would pass no matter what was logged. Attaching
    caplog's own handler to each pipeline logger sidesteps propagation entirely.
    """
    import logging

    names = ["cassa_photometry", "cassa_calibrate", "cassa_integrate",
             "cassa_diagnose", "cassa_run"]
    loggers = [logging.getLogger(n) for n in names]
    for logger in loggers:
        logger.addHandler(caplog.handler)
        logger.setLevel(logging.DEBUG)
    caplog.set_level(logging.DEBUG)
    try:
        yield caplog
    finally:
        for logger in loggers:
            logger.removeHandler(caplog.handler)
