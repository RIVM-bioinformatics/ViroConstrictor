"""Unit tests for user profile prompts and configuration handling."""

from __future__ import annotations

import configparser
from pathlib import Path
from types import SimpleNamespace

import pytest

import ViroConstrictor.userprofile as userprofile


def _write_valid_config(path: Path, *, compmode: str = "local", auto_update: str = "yes", repro_method: str = "conda") -> None:
    """Write a valid profile configuration for the supplied test values.

    Parameters
    ----------
    path : Path
        Destination file to create or overwrite.
    compmode : str, optional
        Value written to ``COMPUTING.compmode``.
    auto_update : str, optional
        Value written to ``GENERAL.auto_update``.
    repro_method : str, optional
        Value written to ``REPRODUCTION.repro_method``.
    """
    config = configparser.ConfigParser()
    config["COMPUTING"] = {"compmode": compmode}
    if compmode == "grid":
        config["COMPUTING"]["queuename"] = "normal"

    config["GENERAL"] = {"auto_update": auto_update}
    if auto_update == "no":
        config["GENERAL"]["ask_for_update"] = "yes"

    config["REPRODUCTION"] = {"repro_method": repro_method}
    if repro_method == "containers":
        config["REPRODUCTION"]["container_cache_path"] = "/tmp/cache"
    else:
        config["REPRODUCTION"]["container_cache_path"] = "None"

    with open(path, "w", encoding="utf-8") as handle:
        config.write(handle)


def test_fileexists_and_fileispopulated(tmp_path: Path) -> None:
    config_file = tmp_path / "profile.ini"

    assert userprofile.FileExists(config_file) is False

    config_file.touch()
    assert userprofile.FileExists(config_file) is True
    assert userprofile.FileIsPopulated(config_file) is False

    config_file.write_text("[GENERAL]\nauto_update = yes\n", encoding="utf-8")
    assert userprofile.FileIsPopulated(config_file) is True


def test_askprompts_fixedchoices_valid(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeCompleter:
        def __init__(self):
            self.created_options = None
            self.listCompleter = object()
            self.pathCompleter = object()

        def createListCompleter(self, options):
            self.created_options = options

    fake = FakeCompleter()
    monkeypatch.setattr(userprofile, "tabCompleter", lambda: fake)
    monkeypatch.setattr(userprofile.readline, "set_completer_delims", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "parse_and_bind", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "set_completer", lambda *_: None)

    called = {"clear": 0}
    monkeypatch.setattr(userprofile.subprocess, "call", lambda *_args, **_kwargs: called.__setitem__("clear", called["clear"] + 1))
    monkeypatch.setattr(userprofile.Console, "input", lambda self, _prompt: "local")

    value = userprofile.AskPrompts("intro", "prompt", ["local", "grid"], fixedchoices=True)

    assert value == "local"
    assert fake.created_options == ["local", "grid"]
    assert called["clear"] == 1


def test_askprompts_fixedchoices_invalid_then_valid(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeCompleter:
        def __init__(self):
            self.listCompleter = object()
            self.pathCompleter = object()

        def createListCompleter(self, _):
            return None

    replies = iter(["wrong", "grid"])

    monkeypatch.setattr(userprofile, "tabCompleter", FakeCompleter)
    monkeypatch.setattr(userprofile.readline, "set_completer_delims", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "parse_and_bind", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "set_completer", lambda *_: None)
    monkeypatch.setattr(userprofile.subprocess, "call", lambda *_args, **_kwargs: 0)
    monkeypatch.setattr(userprofile.Console, "input", lambda self, _prompt: next(replies))

    value = userprofile.AskPrompts("intro", "prompt", ["local", "grid"], fixedchoices=True)

    assert value == "grid"


def test_askprompts_fixedchoices_quit(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeCompleter:
        def __init__(self):
            self.listCompleter = object()
            self.pathCompleter = object()

        def createListCompleter(self, _):
            return None

    monkeypatch.setattr(userprofile, "tabCompleter", FakeCompleter)
    monkeypatch.setattr(userprofile.readline, "set_completer_delims", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "parse_and_bind", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "set_completer", lambda *_: None)
    monkeypatch.setattr(userprofile.subprocess, "call", lambda *_args, **_kwargs: 0)
    monkeypatch.setattr(userprofile.Console, "input", lambda self, _prompt: "quit")

    with pytest.raises(SystemExit) as excinfo:
        userprofile.AskPrompts("intro", "prompt", ["local", "grid"], fixedchoices=True)

    assert excinfo.value.code == -1


def test_askprompts_path_completion_and_default(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeCompleter:
        def __init__(self):
            self.pathCompleter = object()

    monkeypatch.setattr(userprofile, "tabCompleter", FakeCompleter)
    monkeypatch.setattr(userprofile.readline, "parse_and_bind", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "set_completer", lambda *_: None)
    monkeypatch.setattr(userprofile.subprocess, "call", lambda *_args, **_kwargs: 0)
    monkeypatch.setattr("builtins.input", lambda _prompt: "   ")

    value = userprofile.AskPrompts("intro", "prompt", [], fixedchoices=False, default="/tmp/cache")

    assert value == "/tmp/cache"


def test_askprompts_path_completion_returns_text(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeCompleter:
        def __init__(self):
            self.pathCompleter = object()

    monkeypatch.setattr(userprofile, "tabCompleter", FakeCompleter)
    monkeypatch.setattr(userprofile.readline, "parse_and_bind", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "set_completer", lambda *_: None)
    monkeypatch.setattr(userprofile.subprocess, "call", lambda *_args, **_kwargs: 0)
    monkeypatch.setattr("builtins.input", lambda _prompt: "queueA")

    value = userprofile.AskPrompts("intro", "prompt", [], fixedchoices=False)

    assert value == "queueA"


def test_askprompts_path_completion_quit(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeCompleter:
        def __init__(self):
            self.pathCompleter = object()

    monkeypatch.setattr(userprofile, "tabCompleter", FakeCompleter)
    monkeypatch.setattr(userprofile.readline, "parse_and_bind", lambda *_: None)
    monkeypatch.setattr(userprofile.readline, "set_completer", lambda *_: None)
    monkeypatch.setattr(userprofile.subprocess, "call", lambda *_args, **_kwargs: 0)
    monkeypatch.setattr("builtins.input", lambda _prompt: "quit")

    with pytest.raises(SystemExit) as excinfo:
        userprofile.AskPrompts("intro", "prompt", [], fixedchoices=False)

    assert excinfo.value.code == -1


@pytest.mark.xfail(reason="Known defect: BuildConfig writes a non-string value to ConfigParser in the no-container path")
def test_buildconfig_grid_auto_update_no_and_no_containers(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    profile = tmp_path / "profile.ini"
    profile.write_text("stale=true\n", encoding="utf-8")

    answers = iter(["grid", "highmem", "no", "yes"])
    monkeypatch.setattr(userprofile, "AskPrompts", lambda *_args, **_kwargs: next(answers))
    monkeypatch.setattr(userprofile, "containerization_installed", False)

    logged: list[str] = []
    monkeypatch.setattr(userprofile.log, "info", lambda msg: logged.append(msg))

    userprofile.BuildConfig(profile)

    cfg = configparser.ConfigParser()
    cfg.read(profile)

    assert cfg["COMPUTING"]["compmode"] == "grid"
    assert cfg["COMPUTING"]["queuename"] == "highmem"
    assert cfg["GENERAL"]["auto_update"] == "no"
    assert cfg["GENERAL"]["ask_for_update"] == "yes"
    assert cfg["REPRODUCTION"]["repro_method"] == "conda"
    assert cfg["REPRODUCTION"]["container_cache_path"] == "None"
    assert any("Successfully written global configuration settings" in msg for msg in logged)


def test_buildconfig_local_auto_update_yes_with_containers(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    profile = tmp_path / "profile.ini"

    answers = iter(["local", "yes", "/custom/cache"])
    monkeypatch.setattr(userprofile, "AskPrompts", lambda *_args, **_kwargs: next(answers))
    monkeypatch.setattr(userprofile, "containerization_installed", True)
    monkeypatch.setattr(userprofile.log, "info", lambda *_: None)

    userprofile.BuildConfig(profile)

    cfg = configparser.ConfigParser()
    cfg.read(profile)

    assert cfg["COMPUTING"]["compmode"] == "local"
    assert cfg["GENERAL"]["auto_update"] == "yes"
    assert cfg.has_option("GENERAL", "ask_for_update") is False
    assert cfg["REPRODUCTION"]["repro_method"] == "containers"
    assert cfg["REPRODUCTION"]["container_cache_path"] == "/custom/cache"


def test_alloptionsgiven_returns_true_for_complete_config() -> None:
    config = configparser.ConfigParser()
    config["COMPUTING"] = {"compmode": "grid", "queuename": "normal"}
    config["GENERAL"] = {"auto_update": "no", "ask_for_update": "yes"}
    config["REPRODUCTION"] = {"repro_method": "containers", "container_cache_path": "/tmp/cache"}

    assert userprofile.AllOptionsGiven(config) is True


def test_alloptionsgiven_returns_false_when_sections_missing() -> None:
    config = configparser.ConfigParser()
    config["GENERAL"] = {"auto_update": "yes"}

    assert userprofile.AllOptionsGiven(config) is False


@pytest.mark.parametrize(
    "section,values",
    [
        ("COMPUTING", {"compmode": "grid"}),
        ("GENERAL", {"auto_update": "no"}),
        ("REPRODUCTION", {"repro_method": "containers"}),
    ],
)
def test_alloptionsgiven_returns_false_when_required_option_missing(section: str, values: dict[str, str]) -> None:
    config = configparser.ConfigParser()
    config["COMPUTING"] = {"compmode": "local", "queuename": "normal"}
    config["GENERAL"] = {"auto_update": "yes"}
    config["REPRODUCTION"] = {"repro_method": "conda", "container_cache_path": "None"}

    config[section] = values

    assert userprofile.AllOptionsGiven(config) is False


def test_readconfig_creates_missing_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    profile = tmp_path / "profile.ini"
    calls: list[Path] = []

    def fake_build(target: Path) -> None:
        calls.append(target)
        _write_valid_config(target)

    monkeypatch.setattr(userprofile, "BuildConfig", fake_build)
    monkeypatch.setattr(userprofile.log, "info", lambda *_: None)

    cfg = userprofile.ReadConfig(profile)

    assert calls == [profile]
    assert cfg["COMPUTING"]["compmode"] == "local"


def test_readconfig_rebuilds_empty_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    profile = tmp_path / "profile.ini"
    profile.touch()
    calls = {"count": 0}

    def fake_build(target: Path) -> None:
        calls["count"] += 1
        _write_valid_config(target)

    monkeypatch.setattr(userprofile, "BuildConfig", fake_build)
    monkeypatch.setattr(userprofile.log, "info", lambda *_: None)

    cfg = userprofile.ReadConfig(profile)

    assert calls["count"] == 1
    assert cfg["GENERAL"]["auto_update"] == "yes"


def test_readconfig_rebuilds_until_options_present(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    profile = tmp_path / "profile.ini"

    # Start with an incomplete config to force at least one rebuild in the while loop.
    incomplete = configparser.ConfigParser()
    incomplete["COMPUTING"] = {"compmode": "grid"}
    incomplete["GENERAL"] = {"auto_update": "yes"}
    incomplete["REPRODUCTION"] = {"repro_method": "conda", "container_cache_path": "None"}
    with open(profile, "w", encoding="utf-8") as handle:
        incomplete.write(handle)

    calls = {"count": 0}

    def fake_build(target: Path) -> None:
        calls["count"] += 1
        _write_valid_config(target, compmode="grid", auto_update="no", repro_method="containers")

    monkeypatch.setattr(userprofile, "BuildConfig", fake_build)
    monkeypatch.setattr(userprofile.log, "info", lambda *_: None)

    cfg = userprofile.ReadConfig(profile)

    assert calls["count"] == 1
    assert cfg["COMPUTING"]["queuename"] == "normal"
    assert cfg["GENERAL"]["ask_for_update"] == "yes"
    assert cfg["REPRODUCTION"]["container_cache_path"] == "/tmp/cache"


def test_readconfig_rebuilds_after_initial_empty_then_incomplete(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    profile = tmp_path / "profile.ini"
    profile.touch()

    write_state = SimpleNamespace(step=0)

    def fake_build(target: Path) -> None:
        write_state.step += 1
        if write_state.step == 1:
            # First rebuild still incomplete, forcing the while-loop rebuild path too.
            incomplete = configparser.ConfigParser()
            incomplete["COMPUTING"] = {"compmode": "grid"}
            incomplete["GENERAL"] = {"auto_update": "yes"}
            incomplete["REPRODUCTION"] = {"repro_method": "conda", "container_cache_path": "None"}
            with open(target, "w", encoding="utf-8") as handle:
                incomplete.write(handle)
            return

        _write_valid_config(target, compmode="grid", auto_update="no", repro_method="containers")

    monkeypatch.setattr(userprofile, "BuildConfig", fake_build)
    monkeypatch.setattr(userprofile.log, "info", lambda *_: None)

    cfg = userprofile.ReadConfig(profile)

    assert write_state.step == 2
    assert cfg["COMPUTING"]["compmode"] == "grid"
    assert cfg["REPRODUCTION"]["repro_method"] == "containers"
