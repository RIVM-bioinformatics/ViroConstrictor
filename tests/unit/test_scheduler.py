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

        # Test with no scheduler in config
        config.remove_option("COMPUTING", "scheduler")
        assert Scheduler._scheduler_from_config(config) is None

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
        assert (
            Scheduler.determine_scheduler(None, config, dryrun_arg=False)
            is Scheduler.LOCAL
        )

        # Test with scheduler from argument
        assert (
            Scheduler.determine_scheduler("local", config, dryrun_arg=False)
            == Scheduler.LOCAL
        )

        # Test with scheduler from config
        config.set("COMPUTING", "scheduler", "SLURM")
        assert (
            Scheduler.determine_scheduler(None, config, dryrun_arg=False)
            == Scheduler.SLURM
        )

        # Test with environment variable
        with mock.patch.dict("os.environ", {"SLURM_JOB_ID": "123"}):
            assert (
                Scheduler.determine_scheduler(None, config, dryrun_arg=False)
                == Scheduler.SLURM
            )

        # Test with DRMAA
        sys.modules["drmaa"] = fake_drmaa
        assert (
            Scheduler.determine_scheduler(None, config, dryrun_arg=False)
            == Scheduler.SLURM
        )
