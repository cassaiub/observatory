import pytest

from cassa_photometry import cli


@pytest.mark.parametrize("func", ["calibrate", "integrate", "photometry", "verify", "run_all"])
def test_cli_help_exits_cleanly(func, monkeypatch):
    prog = {"calibrate": "cassa-calibrate", "integrate": "cassa-integrate",
            "photometry": "cassa-photometry", "verify": "cassa-verify",
            "run_all": "cassa-run"}[func]
    monkeypatch.setattr("sys.argv", [prog, "--help"])
    with pytest.raises(SystemExit) as exc:
        getattr(cli, func)()
    assert exc.value.code == 0
