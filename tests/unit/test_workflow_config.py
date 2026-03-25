"""Test workflow configuration assembly and resource selection.

These tests verify scheduler-aware configuration, thread and memory limits,
container setup, and resource defaults for core workflow stages.
"""

from argparse import Namespace
from configparser import ConfigParser
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Hashable, cast

import pytest
import yaml
from snakemake_interface_executor_plugins.settings import ExecMode

import ViroConstrictor.workflow_config as workflow_config
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.scheduler import Scheduler


def _base_user_config(*, compmode: str = "local", repro_method: str = "conda", queue: str | None = None) -> ConfigParser:
    """Build a minimal user config for workflow configuration tests.

    Parameters
    ----------
    compmode : str, optional
        Computing mode written to the COMPUTING section.
    repro_method : str, optional
        Reproduction method written to the REPRODUCTION section.
    queue : str | None, optional
        Optional scheduler queue or partition name.

    Returns
    -------
    ConfigParser
        Config parser containing the sections required by WorkflowConfig.
    """
    config = ConfigParser()
    config.add_section("COMPUTING")
    config.set("COMPUTING", "compmode", compmode)
    if queue is not None:
        config.set("COMPUTING", "queuename", queue)

    config.add_section("REPRODUCTION")
    config.set("REPRODUCTION", "repro_method", repro_method)
    config.set("REPRODUCTION", "container_cache_path", "/tmp/cache")
    return config


def _parsed_inputs(tmp_path: Path, *, compmode: str = "local", repro_method: str = "conda", dryrun: bool = False) -> SimpleNamespace:
    """Build parsed input state for WorkflowConfig tests.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory used as the workflow workdir.
    compmode : str, optional
        Computing mode passed through the user config.
    repro_method : str, optional
        Reproduction method passed through the user config.
    dryrun : bool, optional
        Flag value used to exercise dry-run execution paths.

    Returns
    -------
    SimpleNamespace
        Parsed input object with the attributes accessed by WorkflowConfig.
    """
    return SimpleNamespace(
        flags=Namespace(
            dryrun=dryrun,
            platform="illumina",
            unidirectional=True,
            amplicon_type="amplicon",
            threads=8,
        ),
        user_config=_base_user_config(compmode=compmode, repro_method=repro_method),
        match_ref_snakefile="match_ref/workflow.smk",
        snakefile="main/workflow.smk",
        samples_dict={
            "S1": {"INPUTFILE": "single.fastq.gz", "PRESET": "DEFAULT"},
        },
        scheduler=Scheduler.LOCAL,
        workdir=str(tmp_path),
        logfile=str(tmp_path / "viro.log"),
    )


@pytest.fixture
def patch_snakemake_settings(monkeypatch: pytest.MonkeyPatch):
    """Patch Snakemake settings constructors with lightweight stand-ins."""

    def make_ctor(kind: str):
        def ctor(*args, **kwargs):
            return SimpleNamespace(_kind=kind, args=args, **kwargs)

        return ctor

    for name in [
        "ConfigSettings",
        "DAGSettings",
        "DeploymentSettings",
        "ExecutionSettings",
        "OutputSettings",
        "RemoteExecutionSettings",
        "ResourceSettings",
        "SchedulingSettings",
        "StorageSettings",
        "WorkflowSettings",
    ]:
        monkeypatch.setattr(workflow_config, name, make_ctor(name))


def test_correct_unidirectional_flag_returns_true_for_single_end_inputfile() -> None:
    samples: dict[Hashable, Any] = {"S1": {"INPUTFILE": "R1.fastq.gz"}}
    flags = Namespace(platform="illumina", unidirectional=True)
    assert workflow_config.correct_unidirectional_flag(samples, flags) is True


def test_correct_unidirectional_flag_forces_false_when_paired_and_flag_true() -> None:
    samples: dict[Hashable, Any] = {"S1": {"R1": "R1.fastq.gz", "R2": "R2.fastq.gz"}}
    flags = Namespace(platform="illumina", unidirectional=True)
    assert workflow_config.correct_unidirectional_flag(samples, flags) is False


def test_correct_unidirectional_flag_forces_true_when_single_end_and_flag_false() -> None:
    samples: dict[Hashable, Any] = {"S1": {"INPUTFILE": "single.fastq.gz"}}
    flags = Namespace(platform="illumina", unidirectional=False)
    assert workflow_config.correct_unidirectional_flag(samples, flags) is True


def test_correct_unidirectional_flag_returns_false_when_paired_and_flag_false() -> None:
    samples: dict[Hashable, Any] = {"S1": {"R1": "R1.fastq.gz", "R2": "R2.fastq.gz"}}
    flags = Namespace(platform="illumina", unidirectional=False)
    assert workflow_config.correct_unidirectional_flag(samples, flags) is False


def test_correct_unidirectional_flag_defaults_true_for_other_platform() -> None:
    samples: dict[Hashable, Any] = {"S1": {"R1": "x", "R2": "y"}}
    flags = Namespace(platform="nanopore", unidirectional=False)
    assert workflow_config.correct_unidirectional_flag(samples, flags) is True


def test_max_threads_per_type_local_mode_caps_and_minimum() -> None:
    parser = SimpleNamespace(flags=Namespace(threads=2))
    config = _base_user_config(compmode="local")

    result = workflow_config.MaxThreadsPerType(cast(CLIparser, parser), config)

    assert result.assignment is False
    assert 1 <= result.highcpu <= 12
    assert 1 <= result.midcpu <= 6
    assert result.highcpu >= result.midcpu
    assert result.lowcpu == 1


def test_max_threads_per_type_local_mode_scales_with_more_threads() -> None:
    config = _base_user_config(compmode="local")
    small = workflow_config.MaxThreadsPerType(cast(CLIparser, SimpleNamespace(flags=Namespace(threads=4))), config)
    large = workflow_config.MaxThreadsPerType(cast(CLIparser, SimpleNamespace(flags=Namespace(threads=20))), config)

    # More requested threads should not reduce assigned thread classes.
    assert large.highcpu >= small.highcpu
    assert large.midcpu >= small.midcpu
    # Contractually capped to avoid oversubscription in local mode.
    assert large.highcpu <= 12
    assert large.midcpu <= 6
    assert small.lowcpu == 1
    assert large.lowcpu == 1


def test_max_threads_per_type_grid_mode_assigns_fixed_values() -> None:
    parser = SimpleNamespace(flags=Namespace(threads=1))
    config = _base_user_config(compmode="grid")

    result = workflow_config.MaxThreadsPerType(cast(CLIparser, parser), config)

    assert result.assignment is True
    assert result.highcpu == 12
    assert result.midcpu == 6
    assert result.lowcpu == 2


def test_write_yaml_creates_parent_directory_and_writes_content(tmp_path: Path) -> None:
    output = tmp_path / "nested" / "config.yaml"
    payload = {"a": 1, "b": {"x": "y"}}

    returned = workflow_config._write_yaml(payload, str(output))

    assert returned == str(output)
    assert output.exists()
    assert yaml.safe_load(output.read_text()) == payload


@pytest.mark.parametrize(
    "requested,available,expected",
    [
        (8, 8, 6),
        (16, 8, 6),
        (4, 8, 4),
    ],
)
def test_set_cores_behaviour(monkeypatch: pytest.MonkeyPatch, requested: int, available: int, expected: int) -> None:
    obj = workflow_config.WorkflowConfig.__new__(workflow_config.WorkflowConfig)
    monkeypatch.setattr(workflow_config.multiprocessing, "cpu_count", lambda: available)
    assert obj._set_cores(requested) == expected


@pytest.mark.xfail(reason="Current implementation can return <=0 when available CPU count is very low")
def test_set_cores_never_returns_less_than_one(monkeypatch: pytest.MonkeyPatch) -> None:
    obj = workflow_config.WorkflowConfig.__new__(workflow_config.WorkflowConfig)
    monkeypatch.setattr(workflow_config.multiprocessing, "cpu_count", lambda: 1)
    assert obj._set_cores(1) >= 1


def test_get_max_local_mem_uses_sysconf_with_buffer(monkeypatch: pytest.MonkeyPatch) -> None:
    obj = workflow_config.WorkflowConfig.__new__(workflow_config.WorkflowConfig)

    values = {
        "SC_PAGE_SIZE": 1024,
        "SC_PHYS_PAGES": 3000,
    }
    monkeypatch.setattr(workflow_config.os, "sysconf", lambda key: values[key])

    # (1024*3000)/(1024^2) MB ~= 2.93MB; subtract 2000, rounded to nearest 1000 => -2000
    assert obj._get_max_local_mem() == -2000


def test_add_default_resource_settings_returns_empty_for_local() -> None:
    config = _base_user_config(compmode="local", queue="q")
    result = workflow_config.add_default_resource_settings(Scheduler.LOCAL, config)
    # Local execution should not add scheduler-specific resources.
    assert all("lsf_queue" not in arg and "slurm_partition" not in arg for arg in result.args)


def test_add_default_resource_settings_assigns_lsf(monkeypatch: pytest.MonkeyPatch) -> None:
    config = _base_user_config(compmode="grid", queue="normal")

    called: dict[str, Any] = {}

    def fake_assign(default_resource_setting, queue):
        called["queue"] = queue
        called["resource_type"] = type(default_resource_setting).__name__
        return default_resource_setting

    monkeypatch.setattr(workflow_config, "_assign_resources_lsf", fake_assign)
    workflow_config.add_default_resource_settings(Scheduler.LSF, config)
    assert called["queue"] == "normal"
    assert called["resource_type"] == "DefaultResources"


def test_add_default_resource_settings_assigns_slurm(monkeypatch: pytest.MonkeyPatch) -> None:
    config = _base_user_config(compmode="grid", queue="partitionA")

    called: dict[str, Any] = {}

    def fake_assign(default_resource_setting, queue):
        called["queue"] = queue
        called["resource_type"] = type(default_resource_setting).__name__
        return default_resource_setting

    monkeypatch.setattr(workflow_config, "_assign_resources_slurm", fake_assign)
    workflow_config.add_default_resource_settings(Scheduler.SLURM, config)
    assert called["queue"] == "partitionA"
    assert called["resource_type"] == "DefaultResources"


def test_assign_resources_lsf_errors_without_queue() -> None:
    with pytest.raises(SystemExit) as exc:
        workflow_config._assign_resources_lsf(workflow_config.DefaultResources(), None)
    assert exc.value.code == 1


def test_assign_resources_lsf_sets_queue() -> None:
    resources = workflow_config._assign_resources_lsf(workflow_config.DefaultResources(), "short")
    assert "lsf_queue=short" in resources.args


def test_assign_resources_slurm_errors_without_queue() -> None:
    with pytest.raises(SystemExit) as exc:
        workflow_config._assign_resources_slurm(workflow_config.DefaultResources(), None)
    assert exc.value.code == 1


def test_assign_resources_slurm_sets_partition() -> None:
    resources = workflow_config._assign_resources_slurm(workflow_config.DefaultResources(), "gpu")
    assert "slurm_partition=gpu" in resources.args


def test_workflow_config_invalid_stage_raises(tmp_path: Path) -> None:
    parsed = _parsed_inputs(tmp_path)
    with pytest.raises(ValueError, match="VC_stage"):
        workflow_config.WorkflowConfig(parsed_inputs=cast(CLIparser, parsed), vc_stage="INVALID")


def test_workflow_config_main_dryrun_sets_subprocess_execmode(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    patch_snakemake_settings,
) -> None:
    parsed = _parsed_inputs(tmp_path, compmode="local", repro_method="conda", dryrun=True)

    monkeypatch.setattr(workflow_config, "construct_container_bind_args", lambda samples: "bind-args")
    monkeypatch.setattr(workflow_config, "download_containers", lambda prefix, dryrun: 0)
    monkeypatch.setattr(workflow_config.WorkflowConfig, "_get_max_local_mem", lambda self: 64000)
    monkeypatch.setattr(workflow_config.WorkflowConfig, "_set_cores", lambda self, cores: 6)

    cfg = workflow_config.WorkflowConfig(parsed_inputs=cast(CLIparser, parsed), vc_stage="MAIN", outdir_override="od")

    assert cfg.samplesheetfilename == "samples_main"
    assert cfg.workflow_file == "main/workflow.smk"
    assert cfg.output_settings.dryrun is True
    assert cfg.workflow_settings.exec_mode == ExecMode.SUBPROCESS
    assert cfg.resource_settings.cores == 6
    assert cfg.resource_settings.nodes == 1
    assert cfg.snakemake_base_params["outdirOverride"] == "od"
    assert set(cfg.snakemake_base_params["threads"].keys()) == {
        "Alignments",
        "QC",
        "AdapterRemoval",
        "PrimerRemoval",
        "Consensus",
        "Index",
        "Typing",
    }
    assert all(value >= 1 for value in cfg.snakemake_base_params["threads"].values())
    assert cfg.snakemake_base_params["max_local_mem"] == 64000
    assert cfg.snakemake_base_params["computing_execution"] == "local"
    assert cfg.snakemake_base_params["sample_sheet"].endswith("samples_main.yaml")
    assert Path(tmp_path / "config.yaml").exists()
    assert Path(tmp_path / "samples_main.yaml").exists()


def test_workflow_config_mr_grid_container_download_success(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    patch_snakemake_settings,
) -> None:
    parsed = _parsed_inputs(tmp_path, compmode="grid", repro_method="containers", dryrun=False)
    parsed.user_config.set("COMPUTING", "queuename", "queueX")
    parsed.scheduler = Scheduler.SLURM

    called = {}

    monkeypatch.setattr(workflow_config, "construct_container_bind_args", lambda samples: "bind-here")

    def fake_download(prefix, dryrun):
        called["prefix"] = prefix
        called["dryrun"] = dryrun
        return 0

    monkeypatch.setattr(workflow_config, "download_containers", fake_download)
    monkeypatch.setattr(workflow_config.WorkflowConfig, "_get_max_local_mem", lambda self: 50000)

    cfg = workflow_config.WorkflowConfig(parsed_inputs=cast(CLIparser, parsed), vc_stage="MR")

    assert cfg.samplesheetfilename == "samples_mr"
    assert cfg.workflow_file == "match_ref/workflow.smk"
    assert cfg.workflow_settings.exec_mode == ExecMode.DEFAULT
    assert cfg.resource_settings.cores > 1
    assert cfg.resource_settings.nodes > 1
    assert cfg.snakemake_base_params["computing_execution"] == "grid"
    assert all(value >= 1 for value in cfg.snakemake_base_params["threads"].values())
    assert cfg.deployment_settings.apptainer_args == "bind-here"
    assert cfg.deployment_settings.apptainer_prefix == Path("/tmp/cache")
    assert str(called["prefix"]) == "/tmp/cache"
    assert called["dryrun"] is False


def test_workflow_config_container_download_failure_exits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    patch_snakemake_settings,
) -> None:
    parsed = _parsed_inputs(tmp_path, compmode="grid", repro_method="containers", dryrun=False)
    parsed.scheduler = Scheduler.LSF
    parsed.user_config.set("COMPUTING", "queuename", "normal")

    monkeypatch.setattr(workflow_config, "construct_container_bind_args", lambda samples: "bind")
    monkeypatch.setattr(workflow_config, "download_containers", lambda prefix, dryrun: 1)
    monkeypatch.setattr(workflow_config.WorkflowConfig, "_get_max_local_mem", lambda self: 50000)

    with pytest.raises(SystemExit) as exc:
        workflow_config.WorkflowConfig(parsed_inputs=cast(CLIparser, parsed), vc_stage="MAIN")

    assert exc.value.code == 1
