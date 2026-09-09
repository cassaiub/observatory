"""How much of the machine a run may use, and saying so before it starts.

Two things are being pinned. First, that a configured budget is *honoured* --
including the interaction that matters most in practice: a stated budget also
suppresses the interactive CPU prompt, because being asked again would make the
setting useless in exactly the places it is wanted (a notebook, a batch job, a
shared node). Second, that the estimate is available as data, so a notebook can
show a participant what the run needs and what fraction of their machine that
is, rather than reimplementing the arithmetic and drifting from the log line.
"""

from unittest import mock

import pytest

from cassa_photometry.config import load_config
from cassa_photometry.logging_utils import get_logger
from cassa_photometry.phase2_integration.hardware import HardwareManager


def _manager(**budget):
    config = load_config()
    for key, value in budget.items():
        setattr(config.phase2, key, value)
    return HardwareManager(logger=get_logger("t"), config=config)


# --- the budget is honoured ---------------------------------------------------

def test_cpu_fraction_sets_the_core_count():
    hw = _manager(cpu_fraction=0.5)
    cores = hw.allocate_resources()
    assert cores == max(1, round(hw.total_cores * 0.5))


def test_max_cores_caps_whatever_the_fraction_asked_for():
    hw = _manager(cpu_fraction=1.0, max_cores=2)
    assert hw.allocate_resources() == 2


def test_a_configured_budget_does_not_prompt():
    """The whole point of the setting. If this regressed, a notebook cell would
    hang on input() and a batch job would take the default anyway."""
    hw = _manager(cpu_fraction=0.25)
    with mock.patch.object(hw, "_interactive", return_value=True), \
         mock.patch("builtins.input", side_effect=AssertionError("prompted anyway")):
        hw.allocate_resources()


def test_without_a_budget_a_non_interactive_run_takes_the_default():
    hw = _manager()
    with mock.patch.object(hw, "_interactive", return_value=False):
        assert hw.allocate_resources() == max(1, int(hw.total_cores * 0.5))


@pytest.mark.parametrize("fraction", [0.0, -1.0, 5.0])
def test_a_nonsense_fraction_still_yields_a_usable_core_count(fraction):
    hw = _manager(cpu_fraction=fraction)
    assert hw.allocate_resources() >= 1


def test_the_allocation_reaches_the_threading_libraries():
    """Cores nothing acts on are decoration."""
    import os

    hw = _manager(cpu_fraction=1.0, max_cores=3)
    hw.allocate_resources()
    assert os.environ["OMP_NUM_THREADS"] == "3"
    assert os.environ["MKL_NUM_THREADS"] == "3"


# --- the estimate is reportable -----------------------------------------------

def test_the_estimate_says_what_will_be_used_of_what_is_available():
    hw = _manager(cpu_fraction=0.5)
    hw.allocate_resources()
    e = hw.estimate_resources(num_files=60, max_stack=15, total_size_mb=250.0)

    for key in ("cores_used", "cores_total", "cpu_percent", "peak_ram_gb",
                "ram_total_gb", "ram_available_gb", "ram_percent",
                "runtime_s", "fits_in_budget"):
        assert key in e

    assert e["cores_used"] <= e["cores_total"]
    assert 0 < e["cpu_percent"] <= 100.0
    assert e["peak_ram_gb"] > 0


def test_a_memory_budget_smaller_than_the_estimate_is_flagged():
    """The warning has to arrive before the run, not as an OOM kill during it."""
    hw = _manager(cpu_fraction=0.5, max_memory_gb=0.001)
    hw.allocate_resources()
    e = hw.estimate_resources(num_files=60, max_stack=15, total_size_mb=250.0)
    assert e["fits_in_budget"] is False
    assert "phase2.max_memory_gb" in e["ram_limit_label"]


def test_without_a_memory_budget_the_machine_is_the_limit():
    hw = _manager()
    hw.allocate_resources(assume_yes=True)
    e = hw.estimate_resources(num_files=10, max_stack=5, total_size_mb=50.0)
    assert e["ram_limit_gb"] == pytest.approx(hw.sys_ram_gb)
    assert e["ram_limit_label"] == "system memory"


def test_describe_is_readable_before_and_after_an_estimate():
    hw = _manager()
    assert "cores" in hw.describe()          # no estimate yet
    hw.allocate_resources(assume_yes=True)
    hw.estimate_resources(12, 4, 40.0)
    text = hw.describe()
    assert "This machine" in text and "This run will" in text


# --- the flags reach the config ----------------------------------------------

def test_the_cli_flags_override_the_config_file():
    from types import SimpleNamespace

    from cassa_photometry.cli import _apply_resource_args

    config = load_config()
    config.phase2.cpu_fraction = 0.9
    args = SimpleNamespace(cores=4, cpu_fraction=0.25, max_memory_gb=8.0)

    _apply_resource_args(config, args)
    assert config.phase2.max_cores == 4
    assert config.phase2.cpu_fraction == 0.25
    assert config.phase2.max_memory_gb == 8.0
