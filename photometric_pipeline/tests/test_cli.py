import pytest

from cassa_photometry import cli

#: Every console entry point, so a new one cannot be added without its help
#: text being exercised. `diagnose` was missing from this list.
_COMMANDS = {
    "calibrate": "cassa-calibrate", "integrate": "cassa-integrate",
    "photometry": "cassa-photometry", "verify": "cassa-verify",
    "diagnose": "cassa-diagnose", "run_all": "cassa-run",
    "doctor": "cassa-doctor", "simulate": "cassa-simulate",
    "index_fetch": "cassa-index-fetch",
}


@pytest.mark.parametrize("func", sorted(_COMMANDS))
def test_cli_help_exits_cleanly(func, monkeypatch):
    prog = _COMMANDS[func]
    monkeypatch.setattr("sys.argv", [prog, "--help"])
    with pytest.raises(SystemExit) as exc:
        getattr(cli, func)()
    assert exc.value.code == 0


def test_every_console_script_has_a_help_test():
    """A command declared in pyproject but never exercised here is one whose
    argument parsing can break unnoticed."""
    import re

    declared = set(re.findall(r'^(cassa[\w-]*) = "cassa_photometry\.cli:',
                              open("pyproject.toml").read(), flags=re.MULTILINE))
    declared.discard("cassa")  # the umbrella, covered separately below
    assert declared == set(_COMMANDS.values())


def test_the_umbrella_lists_every_subcommand():
    assert cli.main([]) == 0
    for name in cli.COMMANDS:
        assert name


def test_the_umbrella_rejects_an_unknown_subcommand(capsys):
    assert cli.main(["nope"]) == 2
    assert "unknown command" in capsys.readouterr().err


# --- Step-plan flags (--skip / --only / --show-plan / --from / --to) ----------

def _run(monkeypatch, capsys, argv, func="calibrate"):
    """Invoke an entry point with argv, returning (exit code, stdout)."""
    monkeypatch.setattr("sys.argv", argv)
    with pytest.raises(SystemExit) as exc:
        getattr(cli, func)()
    return exc.value.code, capsys.readouterr().out


def test_show_plan_prints_the_order_and_exits(monkeypatch, capsys):
    code, out = _run(monkeypatch, capsys,
                     ["cassa-calibrate", "-i", "x", "-o", "y", "--show-plan"])
    assert code == 0
    assert "1. linearity" in out
    assert "bias" in out and "flat" in out


def test_skip_removes_a_step_from_the_plan(monkeypatch, capsys):
    code, out = _run(monkeypatch, capsys,
                     ["cassa-calibrate", "-i", "x", "-o", "y",
                      "--skip", "flat", "--show-plan"])
    assert code == 0
    assert "flat   [off]" in out


def test_only_keeps_just_the_named_steps(monkeypatch, capsys):
    code, out = _run(monkeypatch, capsys,
                     ["cassa-calibrate", "-i", "x", "-o", "y",
                      "--only", "bias", "--show-plan"])
    assert code == 0
    assert "1. bias" in out
    assert "dark   [off]" in out


def test_an_unknown_step_name_is_rejected(monkeypatch, capsys):
    code, _ = _run(monkeypatch, capsys,
                   ["cassa-calibrate", "-i", "x", "-o", "y", "--skip", "flatt"])
    assert code != 0


def test_an_ambiguous_bare_step_name_is_rejected(monkeypatch, capsys):
    """`measure_fwhm` exists in both phase 1 and phase 2, so on cassa-run a bare
    name must be an error rather than a guess at which was meant."""
    code, _ = _run(monkeypatch, capsys,
                   ["cassa-run", "-i", "x", "-o", "y", "--skip", "measure_fwhm"],
                   func="run_all")
    assert code != 0


def test_a_qualified_step_name_resolves(monkeypatch, capsys):
    code, out = _run(monkeypatch, capsys,
                     ["cassa-run", "-i", "x", "-o", "y",
                      "--skip", "phase2.measure_fwhm", "--show-plan"],
                     func="run_all")
    assert code == 0
    assert "measure_fwhm   [off]" in out
    # Phase 1's identically-named step is untouched.
    assert "6. measure_fwhm" in out or "8. measure_fwhm" in out


def test_skipping_a_step_another_needs_is_still_allowed(monkeypatch, capsys):
    """Excluding is always legal: if nothing provides the token, nothing
    requires it."""
    code, out = _run(monkeypatch, capsys,
                     ["cassa-calibrate", "-i", "x", "-o", "y",
                      "--skip", "bias", "--show-plan"])
    assert code == 0
    assert "bias   [off]" in out
    assert "dark" in out


def test_from_after_to_is_rejected(monkeypatch, capsys):
    code, _ = _run(monkeypatch, capsys,
                   ["cassa-run", "-i", "x", "-o", "y", "--from", "3", "--to", "1"],
                   func="run_all")
    assert code != 0


def test_run_show_plan_covers_all_three_phases(monkeypatch, capsys):
    code, out = _run(monkeypatch, capsys,
                     ["cassa-run", "-i", "x", "-o", "y", "--show-plan"],
                     func="run_all")
    assert code == 0
    for phase in ("phase1:", "phase2:", "phase3:"):
        assert phase in out
