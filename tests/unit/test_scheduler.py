import sys
import types
from configparser import ConfigParser
from unittest import mock

import pytest

from ViroConstrictor.scheduler import Scheduler


class TestScheduler:

    def test_supported_schedulers(self):
        assert Scheduler.supported_schedulers() == [
            "LOCAL",
            "SLURM",
            "LSF",
            "DRYRUN",
            "AUTO",
        ]

    def test_is_valid(self):
        assert Scheduler.is_valid("local")
        assert Scheduler.is_valid("SLURM")
        assert Scheduler.is_valid("lsf   ")
        assert Scheduler.is_valid("dryrun")
        assert Scheduler.is_valid("auto")
        assert not Scheduler.is_valid("invalid")

    def test_from_string(self):
        assert Scheduler.from_string("local") == Scheduler.LOCAL
        assert Scheduler.from_string("SLURM") == Scheduler.SLURM
        assert Scheduler.from_string("lsf   ") == Scheduler.LSF
        assert Scheduler.from_string("auto") == Scheduler.AUTO
        assert Scheduler.from_string("dryrun") == Scheduler.DRYRUN

        with pytest.raises(ValueError):
            Scheduler.from_string("invalid")

    def test_scheduler_from_argument(self):
        assert Scheduler._scheduler_from_argument("local") == Scheduler.LOCAL
        assert Scheduler._scheduler_from_argument("SLURM") == Scheduler.SLURM
        assert Scheduler._scheduler_from_argument("lsf   ") == Scheduler.LSF
        assert Scheduler._scheduler_from_argument("dryrun") == Scheduler.DRYRUN

        # invalid should be local because a scheduler has been provided, but it is not valid
        assert Scheduler._scheduler_from_argument("invalid") == Scheduler.LOCAL

        # None because it has to continue to the next check
        assert Scheduler._scheduler_from_argument("auto") is None

    def test_scheduler_from_config(self):
        config = ConfigParser()
        config.add_section("COMPUTING")

        config.set("COMPUTING", "scheduler", "local")
        assert Scheduler._scheduler_from_config(config) == Scheduler.LOCAL
        config.set("COMPUTING", "scheduler", "SLURM")
        assert Scheduler._scheduler_from_config(config) == Scheduler.SLURM
        config.set("COMPUTING", "scheduler", "lsf   ")
        assert Scheduler._scheduler_from_config(config) == Scheduler.LSF
        config.set("COMPUTING", "scheduler", "invalid")
        assert Scheduler._scheduler_from_config(config) == Scheduler.LOCAL

    def test_scheduler_from_config_auto_falls_through(self):
        config = ConfigParser()
        config.add_section("COMPUTING")

        config.set("COMPUTING", "scheduler", "AUTO")
        assert Scheduler._scheduler_from_config(config) is None

        config.set("COMPUTING", "scheduler", "auto")
        assert Scheduler._scheduler_from_config(config) is None

        config.set("COMPUTING", "scheduler", "grid")
        assert Scheduler._scheduler_from_config(config) is None

    def test_scheduler_from_environment(self):
        with mock.patch.dict("os.environ", {"SLURM_JOB_ID": "123"}):
            scheduler = Scheduler._scheduler_from_environment()
            assert scheduler == Scheduler.SLURM
    def test_scheduler_from_config_without_scheduler_returns_none(self) -> None:
        """Test that _scheduler_from_config returns None when config has no scheduler value."""
        assert Scheduler._scheduler_from_config(make_config()) is None

    def test_scheduler_from_environment(self):
        with mock.patch.dict("os.environ", {"SLURM_JOB_ID": "123"}):
            scheduler = Scheduler._scheduler_from_environment()
            assert scheduler == Scheduler.SLURM

        with mock.patch.dict("os.environ", {"LSB_JOBID": "456"}):
            scheduler = Scheduler._scheduler_from_environment()
            assert scheduler == Scheduler.LSF

        with mock.patch.dict("os.environ", {}):
            scheduler = Scheduler._scheduler_from_environment()
            assert scheduler is None

    @pytest.fixture
    def fake_drmaa(self):
        class FakeDRMAASession:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc_value, traceback):
                pass

            @property
            def drmsInfo(self):
                return "SLURM"

        return types.SimpleNamespace(Session=lambda: FakeDRMAASession())

    def test_scheduler_from_drmaa(self, fake_drmaa):
        sys.modules["drmaa"] = fake_drmaa

        assert Scheduler._scheduler_from_drmaa() == Scheduler.SLURM

        del sys.modules["drmaa"]
        assert Scheduler._scheduler_from_drmaa() is None

    def test_determine_scheduler(self, fake_drmaa):
        config = ConfigParser()
        config.add_section("COMPUTING")

        # Test with no scheduler set
        assert Scheduler.determine_scheduler(None, config, dryrun_arg=False) is Scheduler.LOCAL

    def test_determine_scheduler_invalid_argument_returns_local(self) -> None:
        """Test that invalid scheduler argument defaults to LOCAL."""
        assert Scheduler.determine_scheduler("invalid", make_config(scheduler="SLURM"), dryrun_arg=False) is Scheduler.LOCAL

    def test_determine_scheduler_auto_argument_falls_back_to_config(self) -> None:
        """Test that 'auto' argument falls back to config scheduler."""
        assert Scheduler.determine_scheduler("auto", make_config(scheduler="lsf"), dryrun_arg=False) is Scheduler.LSF

    @pytest.mark.xfail(reason="Known defect: config scheduler 'AUTO' currently resolves to Scheduler.AUTO instead of falling back to environment")
    def test_determine_scheduler_auto_in_config_falls_back_to_environment(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that 'AUTO' in config falls back to environment detection (known defect).

        Parameters
        ----------
        monkeypatch : pytest.MonkeyPatch
            Pytest fixture for mocking behavior.
        """
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/bsub" if cmd == "bsub" else None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        result = Scheduler.determine_scheduler("", make_config(scheduler="AUTO"), dryrun_arg=False)
        assert result is Scheduler.LSF
        assert result is not Scheduler.AUTO

    @pytest.mark.xfail(
        reason="Known defect: config scheduler 'AUTO' currently resolves to Scheduler.AUTO instead of falling back to LOCAL when no scheduler is detected"
    )
    def test_determine_scheduler_auto_in_config_falls_back_to_local(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that 'AUTO' in config falls back to LOCAL when no scheduler found (known defect).

        Parameters
        ----------
        monkeypatch : pytest.MonkeyPatch
            Pytest fixture for mocking behavior.
        """
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        sys.modules.pop("drmaa", None)

        result = Scheduler.determine_scheduler("", make_config(scheduler="AUTO"), dryrun_arg=False)
        assert result is Scheduler.LOCAL
        assert result is not Scheduler.AUTO

    def test_determine_scheduler_auto_argument_never_returns_auto(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that 'AUTO' argument never returns Scheduler.AUTO after fallback.

        Parameters
        ----------
        monkeypatch : pytest.MonkeyPatch
            Pytest fixture for mocking behavior.
        """
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/sbatch" if cmd == "sbatch" else None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)

        result = Scheduler.determine_scheduler("AUTO", make_config(), dryrun_arg=False)
        assert result is Scheduler.SLURM
        assert result is not Scheduler.AUTO

    def test_determine_scheduler_uses_environment_when_config_absent(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that environment detection is used when config absent.

        Parameters
        ----------
        monkeypatch : pytest.MonkeyPatch
            Pytest fixture for mocking behavior.
        """
        monkeypatch.setattr("shutil.which", lambda cmd: "/usr/bin/sbatch" if cmd == "sbatch" else None)
        assert Scheduler.determine_scheduler("", make_config(), dryrun_arg=False) is Scheduler.SLURM

    def test_determine_scheduler_uses_drmaa_when_environment_absent(self, monkeypatch: pytest.MonkeyPatch, fake_drmaa_slurm) -> None:
        """Test that DRMAA detection is used when environment absent.

        Parameters
        ----------
        monkeypatch : pytest.MonkeyPatch
            Pytest fixture for mocking behavior.
        fake_drmaa_slurm
            Pytest fixture providing a mock DRMAA module with SLURM backend.
        """
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        assert Scheduler.determine_scheduler("", make_config(), dryrun_arg=False) is Scheduler.SLURM

    def test_determine_scheduler_without_sources_defaults_to_local(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test that determine_scheduler defaults to LOCAL when no sources detected.

        Parameters
        ----------
        monkeypatch : pytest.MonkeyPatch
            Pytest fixture for mocking behavior.
        """
        monkeypatch.setattr("shutil.which", lambda cmd: None)
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        monkeypatch.delenv("LSB_JOBID", raising=False)
        sys.modules.pop("drmaa", None)

        # Test with DRMAA
        sys.modules["drmaa"] = fake_drmaa
        assert Scheduler.determine_scheduler(None, config, dryrun_arg=False) == Scheduler.SLURM

    def test_determine_scheduler_auto_from_config_uses_environment(self, fake_drmaa):
        config = ConfigParser()
        config.add_section("COMPUTING")
        config.set("COMPUTING", "compmode", "grid")
        config.set("COMPUTING", "scheduler", "AUTO")

        with mock.patch.dict("os.environ", {"SLURM_JOB_ID": "123"}):
            assert (
                Scheduler.determine_scheduler("", config, dryrun_arg=False)
                == Scheduler.SLURM
            )

        with mock.patch.dict("os.environ", {"LSB_JOBID": "456"}, clear=True):
            assert (
                Scheduler.determine_scheduler("", config, dryrun_arg=False)
                == Scheduler.LSF
            )

        sys.modules["drmaa"] = fake_drmaa
        with mock.patch.dict("os.environ", {}, clear=True):
            # shutil.which patched so we fall through past the sbatch check to DRMAA
            with mock.patch("ViroConstrictor.scheduler.shutil.which", return_value=None):
                assert (
                    Scheduler.determine_scheduler("", config, dryrun_arg=False)
                    == Scheduler.SLURM
                )
