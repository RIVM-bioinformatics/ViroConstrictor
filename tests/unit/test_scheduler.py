"""Test scheduler discovery and precedence rules.

These tests cover scheduler parsing from arguments, config files, runtime
environment signals, and DRMAA-backed detection.
"""

import sys
import types
from configparser import ConfigParser

import pytest

from ViroConstrictor.scheduler import Scheduler


def make_config(*, scheduler: str = "", compmode: str = "") -> ConfigParser:
    """Build a minimal computing config for scheduler tests.

    Parameters
    ----------
    scheduler : str, optional
        Scheduler value written to the COMPUTING section.
    compmode : str, optional
        Computing mode written to the COMPUTING section.

    Returns
    -------
    ConfigParser
        Config parser containing the requested computing settings.
    """
    config = ConfigParser()
    config.add_section("COMPUTING")
    if scheduler:
        config.set("COMPUTING", "scheduler", scheduler)
    if compmode:
        config.set("COMPUTING", "compmode", compmode)
    return config


@pytest.fixture
def fake_drmaa_slurm(monkeypatch):
    """Provide a fake DRMAA session that reports SLURM."""

    class FakeDRMAASession:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        @property
        def drmsInfo(self):
            return "slurm"

    module = types.SimpleNamespace(Session=lambda: FakeDRMAASession())
    monkeypatch.setitem(sys.modules, "drmaa", module)
    return module


@pytest.fixture
def fake_drmaa_invalid(monkeypatch):
    """Provide a fake DRMAA session that reports an unsupported backend."""

    class FakeDRMAASession:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        @property
        def drmsInfo(self):
            return "pbs"

    module = types.SimpleNamespace(Session=lambda: FakeDRMAASession())
    monkeypatch.setitem(sys.modules, "drmaa", module)
    return module


@pytest.fixture
def fake_drmaa_runtime_error(monkeypatch):
    """Provide a fake DRMAA session that raises when entered."""

    class FailingSession:
        def __enter__(self):
            raise RuntimeError("no drmaa backend")

        def __exit__(self, exc_type, exc_value, traceback):
            return False

    module = types.SimpleNamespace(Session=lambda: FailingSession())
    monkeypatch.setitem(sys.modules, "drmaa", module)
    return module


class TestScheduler:
    """Group tests for scheduler detection and precedence behavior."""

    def test_supported_schedulers(self):
        assert Scheduler.supported_schedulers() == ["LOCAL", "SLURM", "LSF", "DRYRUN", "AUTO"]

    @pytest.mark.parametrize(
        "scheduler_name,expected",
        [
            ("local", True),
            ("SLURM", True),
            ("lsf   ", True),
            ("dryrun", True),
            ("auto", True),
            ("grid", True),
            ("none", True),
            ("", True),
            ("invalid", False),
        ],
    )
    def test_is_valid(self, scheduler_name, expected):
        assert Scheduler.is_valid(scheduler_name) is expected

    @pytest.mark.parametrize(
        "scheduler_name,expected",
        [
            ("local", Scheduler.LOCAL),
            ("SLURM", Scheduler.SLURM),
            ("lsf   ", Scheduler.LSF),
            ("auto", Scheduler.AUTO),
            ("grid", Scheduler.AUTO),
            ("dryrun", Scheduler.DRYRUN),
        ],
    )
    def test_from_string(self, scheduler_name, expected):
        assert Scheduler.from_string(scheduler_name) is expected

    def test_from_string_invalid_raises_value_error(self):
        with pytest.raises(ValueError, match="Invalid scheduler"):
            Scheduler.from_string("invalid")

    @pytest.mark.parametrize(
        "scheduler_name,expected",
        [
            ("local", Scheduler.LOCAL),
            ("SLURM", Scheduler.SLURM),
            ("lsf   ", Scheduler.LSF),
            ("dryrun", Scheduler.DRYRUN),
            ("", None),
            (None, None),
        ],
    )
    def test_scheduler_from_argument(self, scheduler_name, expected):
        assert Scheduler._scheduler_from_argument(scheduler_name) is expected

    def test_scheduler_from_argument_invalid_returns_local_and_logs_warning(self, capsys):
        scheduler = Scheduler._scheduler_from_argument("invalid")
        assert scheduler is Scheduler.LOCAL
        assert "Invalid scheduler string" in capsys.readouterr().err

    def test_scheduler_from_argument_auto_returns_none(self):
        assert Scheduler._scheduler_from_argument("auto") is None

    def test_scheduler_from_config_compmode_local_takes_precedence(self):
        config = make_config(scheduler="SLURM", compmode="local")
        assert Scheduler._scheduler_from_config(config) is Scheduler.LOCAL

    def test_scheduler_from_config_invalid_returns_local_and_logs_warning(self, capsys):
        config = make_config(scheduler="invalid")
        scheduler = Scheduler._scheduler_from_config(config)
        assert scheduler is Scheduler.LOCAL
        assert "Invalid scheduler in config" in capsys.readouterr().err

    @pytest.mark.parametrize(
        "scheduler_name,expected",
        [("local", Scheduler.LOCAL), ("SLURM", Scheduler.SLURM), ("lsf   ", Scheduler.LSF)],
    )
    def test_scheduler_from_config_valid(self, scheduler_name, expected):
        assert Scheduler._scheduler_from_config(make_config(scheduler=scheduler_name)) is expected

    def test_scheduler_from_config_without_scheduler_returns_none(self):
        assert Scheduler._scheduler_from_config(make_config()) is None

    def test_scheduler_from_environment_uses_sbatch_binary(self, monkeypatch):
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/sbatch" if cmd == "sbatch" else None)
        assert Scheduler._scheduler_from_environment() is Scheduler.SLURM

    def test_scheduler_from_environment_uses_bsub_binary(self, monkeypatch):
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/bsub" if cmd == "bsub" else None)
        assert Scheduler._scheduler_from_environment() is Scheduler.LSF

    def test_scheduler_from_environment_uses_environment_variables(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        monkeypatch.setenv("SLURM_JOB_ID", "123")
        assert Scheduler._scheduler_from_environment() is Scheduler.SLURM

        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.setenv("LSB_JOBID", "456")
        assert Scheduler._scheduler_from_environment() is Scheduler.LSF

    def test_scheduler_from_environment_without_signals_returns_none(self, monkeypatch):
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        assert Scheduler._scheduler_from_environment() is None

    def test_scheduler_from_drmaa_valid(self, fake_drmaa_slurm):
        assert Scheduler._scheduler_from_drmaa() is Scheduler.SLURM

    def test_scheduler_from_drmaa_runtime_error_returns_none(self, fake_drmaa_runtime_error):
        assert Scheduler._scheduler_from_drmaa() is None

    def test_scheduler_from_drmaa_unknown_scheduler_raises_value_error(self, fake_drmaa_invalid):
        with pytest.raises(ValueError, match="Invalid scheduler"):
            Scheduler._scheduler_from_drmaa()

    def test_determine_scheduler_dryrun_has_highest_precedence(self):
        assert Scheduler.determine_scheduler("SLURM", make_config(scheduler="lsf"), dryrun_arg=True) is Scheduler.DRYRUN

    def test_determine_scheduler_uses_scheduler_argument(self):
        assert Scheduler.determine_scheduler("local", make_config(scheduler="SLURM"), dryrun_arg=False) is Scheduler.LOCAL

    def test_determine_scheduler_invalid_argument_returns_local(self):
        assert Scheduler.determine_scheduler("invalid", make_config(scheduler="SLURM"), dryrun_arg=False) is Scheduler.LOCAL

    def test_determine_scheduler_auto_argument_falls_back_to_config(self):
        assert Scheduler.determine_scheduler("auto", make_config(scheduler="lsf"), dryrun_arg=False) is Scheduler.LSF

    @pytest.mark.xfail(reason="Known defect: config scheduler 'AUTO' currently resolves to Scheduler.AUTO instead of falling back to environment")
    def test_determine_scheduler_auto_in_config_falls_back_to_environment(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/bsub" if cmd == "bsub" else None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        result = Scheduler.determine_scheduler("", make_config(scheduler="AUTO"), dryrun_arg=False)
        assert result is Scheduler.LSF
        assert result is not Scheduler.AUTO

    @pytest.mark.xfail(
        reason="Known defect: config scheduler 'AUTO' currently resolves to Scheduler.AUTO instead of falling back to LOCAL when no scheduler is detected"
    )
    def test_determine_scheduler_auto_in_config_falls_back_to_local(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        sys.modules.pop("drmaa", None)

        result = Scheduler.determine_scheduler("", make_config(scheduler="AUTO"), dryrun_arg=False)
        assert result is Scheduler.LOCAL
        assert result is not Scheduler.AUTO

    def test_determine_scheduler_auto_argument_never_returns_auto(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/sbatch" if cmd == "sbatch" else None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)

        result = Scheduler.determine_scheduler("AUTO", make_config(), dryrun_arg=False)
        assert result is Scheduler.SLURM
        assert result is not Scheduler.AUTO

    def test_determine_scheduler_uses_environment_when_config_absent(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/sbatch" if cmd == "sbatch" else None)
        assert Scheduler.determine_scheduler("", make_config(), dryrun_arg=False) is Scheduler.SLURM

    def test_determine_scheduler_uses_drmaa_when_environment_absent(self, monkeypatch, fake_drmaa_slurm):
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        assert Scheduler.determine_scheduler("", make_config(), dryrun_arg=False) is Scheduler.SLURM

    def test_determine_scheduler_without_sources_defaults_to_local(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        sys.modules.pop("drmaa", None)

        assert Scheduler.determine_scheduler("", "", dryrun_arg=False) is Scheduler.LOCAL
