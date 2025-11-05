from typing import Any
from unittest.mock import MagicMock

import pytest


class TestMemoryAllocation:
    """Test suite for memory allocation functions used in Snakemake workflows."""

    @pytest.fixture
    def mock_wildcards(self):
        """Create mock wildcards object as used in Snakemake."""
        wildcards = MagicMock()
        wildcards.sample = "test_sample"
        wildcards.RefID = "test_ref"
        wildcards.Virus = "test_virus"
        return wildcards

    def low_memory_job(self, wildcards: Any, threads: int, attempt: int, config: dict[str, [str | int]]) -> int:
        """Low memory job function for testing."""
        if config["computing_execution"] == "local":
            return min(attempt * threads * 1 * 1000, int(config["max_local_mem"]))
        return attempt * threads * 1 * 1000

    def medium_memory_job(self, wildcards: Any, threads: int, attempt: int, config: dict[str, [str | int]]) -> int:
        """Medium memory job function for testing."""
        if config["computing_execution"] == "local":
            return min(attempt * threads * 2 * 1000, int(config["max_local_mem"]))
        return attempt * threads * 2 * 1000

    def high_memory_job(self, wildcards: Any, threads: int, attempt: int, config: dict[str, [str | int]]) -> int:
        """High memory job function for testing."""
        if config["computing_execution"] == "local":
            return min(attempt * threads * 4 * 1000, int(config["max_local_mem"]))
        return attempt * threads * 4 * 1000

    def test_low_memory_job_local_execution(self, mock_wildcards):
        """Test low_memory_job function with local execution."""
        config = {"computing_execution": "local", "max_local_mem": 8000}

        result = self.low_memory_job(mock_wildcards, threads=2, attempt=1, config=config)
        assert result == 2000  # 1 * 2 * 1 * 1000

        result = self.low_memory_job(mock_wildcards, threads=2, attempt=3, config=config)
        assert result == 6000  # 3 * 2 * 1 * 1000

        # test cap enforcement
        result = self.low_memory_job(mock_wildcards, threads=4, attempt=3, config=config)
        assert result == 8000  # min(3 * 4 * 1 * 1000, 8000) = min(12000, 8000)

    def test_low_memory_job_cluster_execution(self, mock_wildcards):
        """Test low_memory_job function with cluster execution (no memory cap)."""
        config = {"computing_execution": "cluster", "max_local_mem": 8000}

        result = self.low_memory_job(mock_wildcards, threads=4, attempt=3, config=config)
        assert result == 12000  # 3 * 4 * 1 * 1000

        # should not cap at 8000
        result = self.low_memory_job(mock_wildcards, threads=8, attempt=5, config=config)
        assert result == 40000  # 5 * 8 * 1 * 1000

    def test_medium_memory_job_local_execution(self, mock_wildcards):
        """Test medium_memory_job function with local execution."""
        config = {"computing_execution": "local", "max_local_mem": 16000}

        result = self.medium_memory_job(mock_wildcards, threads=2, attempt=1, config=config)
        assert result == 4000  # 1 * 2 * 2 * 1000

        result = self.medium_memory_job(mock_wildcards, threads=2, attempt=3, config=config)
        assert result == 12000  # 3 * 2 * 2 * 1000

        # Test memory cap enforcement
        result = self.medium_memory_job(mock_wildcards, threads=4, attempt=3, config=config)
        assert result == 16000  # min(3 * 4 * 2 * 1000, 16000) = min(24000, 16000)

    def test_medium_memory_job_cluster_execution(self, mock_wildcards):
        """Test medium_memory_job function with cluster execution."""
        config = {"computing_execution": "cluster", "max_local_mem": 16000}

        result = self.medium_memory_job(mock_wildcards, threads=4, attempt=3, config=config)
        assert result == 24000  # 3 * 4 * 2 * 1000 (no cap applied)

    def test_high_memory_job_local_execution(self, mock_wildcards):
        """Test high_memory_job function with local execution."""
        config = {"computing_execution": "local", "max_local_mem": 32000}

        result = self.high_memory_job(mock_wildcards, threads=2, attempt=1, config=config)
        assert result == 8000  # 1 * 2 * 4 * 1000

        result = self.high_memory_job(mock_wildcards, threads=2, attempt=3, config=config)
        assert result == 24000  # 3 * 2 * 4 * 1000

        # Test memory cap enforcement
        result = self.high_memory_job(mock_wildcards, threads=4, attempt=3, config=config)
        assert result == 32000  # min(3 * 4 * 4 * 1000, 32000) = min(48000, 32000)

    def test_high_memory_job_cluster_execution(self, mock_wildcards):
        """Test high_memory_job function with cluster execution."""
        config = {"computing_execution": "cluster", "max_local_mem": 32000}

        # Test basic calculation without cap
        result = self.high_memory_job(mock_wildcards, threads=4, attempt=3, config=config)
        assert result == 48000  # 3 * 4 * 4 * 1000 (no cap applied)

    def test_zero_threads_edge_case(self, mock_wildcards):
        """Test behavior with zero threads."""
        config = {"computing_execution": "local", "max_local_mem": 8000}

        # With 0 threads, all functions should return 0
        low_result = self.low_memory_job(mock_wildcards, 0, 1, config)
        medium_result = self.medium_memory_job(mock_wildcards, 0, 1, config)
        high_result = self.high_memory_job(mock_wildcards, 0, 1, config)

        assert low_result == 0
        assert medium_result == 0
        assert high_result == 0

    def test_function_works_as_snakemake_resource(self, mock_wildcards):
        """Test that functions work when used as Snakemake resource functions."""
        config = {"computing_execution": "local", "max_local_mem": 8000}

        # Simulate how Snakemake calls these functions
        # (with keyword arguments and as resource allocation)

        # Test that functions can be called with keyword arguments
        result = self.low_memory_job(wildcards=mock_wildcards, threads=2, attempt=1, config=config)
        assert result == 2000

        # Test that result can be used as memory allocation (positive integer)
        assert result > 0
        assert isinstance(result, int)

    def test_retry_escalation_behavior(self, mock_wildcards):
        """Test that memory increases with retry attempts."""
        config = {"computing_execution": "cluster", "max_local_mem": 100000}

        threads = 2

        # Test escalation for low memory jobs
        attempt1 = self.low_memory_job(mock_wildcards, threads, 1, config)
        attempt2 = self.low_memory_job(mock_wildcards, threads, 2, config)
        attempt3 = self.low_memory_job(mock_wildcards, threads, 3, config)

        assert attempt2 == 2 * attempt1
        assert attempt3 == 3 * attempt1

        # Test escalation for medium memory jobs
        attempt1_medium = self.medium_memory_job(mock_wildcards, threads, 1, config)
        attempt2_medium = self.medium_memory_job(mock_wildcards, threads, 2, config)
        attempt3_medium = self.medium_memory_job(mock_wildcards, threads, 3, config)
        assert attempt2_medium == 2 * attempt1_medium
        assert attempt3_medium == 3 * attempt1_medium

        # Test escalation for high memory jobs
        attempt1_high = self.high_memory_job(mock_wildcards, threads, 1, config)
        attempt2_high = self.high_memory_job(mock_wildcards, threads, 2, config)
        attempt3_high = self.high_memory_job(mock_wildcards, threads, 3, config)
        assert attempt2_high == 2 * attempt1_high
        assert attempt3_high == 3 * attempt1_high

    def test_snakemake_rule_resource_usage_simulation(self, mock_wildcards):
        """Simulate how these functions would be used in actual Snakemake rules."""

        # Simulate rule resource definition like: resources: mem_mb=high_memory_job
        # This is how the functions are actually called in the workflow

        class MockSnakemakeRule:
            def __init__(self, memory_func):
                self.memory_func = memory_func

            def allocate_memory(self, wildcards, threads, attempt):
                """Simulate Snakemake calling the memory function."""
                return self.memory_func(wildcards, threads, attempt)

        # Create rules with our memory functions (simulating workflow.smk usage)
        low_mem_rule = MockSnakemakeRule(
            lambda wc, threads, attempt: self.low_memory_job(wc, threads, attempt, {"computing_execution": "local", "max_local_mem": 8000})
        )

        high_mem_rule = MockSnakemakeRule(
            lambda wc, threads, attempt: self.high_memory_job(wc, threads, attempt, {"computing_execution": "cluster", "max_local_mem": 32000})
        )

        # Test rule execution
        low_result = low_mem_rule.allocate_memory(mock_wildcards, 2, 1)
        high_result = high_mem_rule.allocate_memory(mock_wildcards, 2, 1)

        assert low_result == 2000
        assert high_result == 8000

        # Test retry behavior
        retry_result = low_mem_rule.allocate_memory(mock_wildcards, 2, 3)
        assert retry_result == 6000  # 3x the initial memory request


class TestRuntimeAllocation:
    """Test suite for runtime allocation functions used in Snakemake workflows."""

    @pytest.fixture
    def mock_wildcards(self):
        """Create mock wildcards object as used in Snakemake."""
        wildcards = MagicMock()
        wildcards.sample = "test_sample"
        wildcards.RefID = "test_ref"
        wildcards.Virus = "test_virus"
        return wildcards

    def low_runtime_job(self, wildcards: Any, attempt: int) -> int:
        return attempt * 2

    def medium_runtime_job(self, wildcards: Any, attempt: int) -> int:
        return attempt * 10

    def high_runtime_job(self, wildcards: Any, attempt: int) -> int:
        return attempt * 30

    def test_low_runtime_job(self, mock_wildcards):
        """Test low_runtime_job function."""
        result = self.low_runtime_job(mock_wildcards, attempt=1)
        assert result == 2  # 1 * 2

        result = self.low_runtime_job(mock_wildcards, attempt=3)
        assert result == 6  # 3 * 2

    def test_medium_runtime_job(self, mock_wildcards):
        """Test medium_runtime_job function."""
        result = self.medium_runtime_job(mock_wildcards, attempt=1)
        assert result == 10  # 1 * 10

        result = self.medium_runtime_job(mock_wildcards, attempt=3)
        assert result == 30  # 3 * 10

    def test_high_runtime_job(self, mock_wildcards):
        """Test high_runtime_job function."""
        result = self.high_runtime_job(mock_wildcards, attempt=1)
        assert result == 30  # 1 * 30

        result = self.high_runtime_job(mock_wildcards, attempt=3)
        assert result == 90  # 3 * 30

    def test_runtime_functions_as_snakemake_resources(self, mock_wildcards):
        """Test that runtime functions work when used as Snakemake resource functions."""
        # Test that functions can be called with keyword arguments
        low_result = self.low_runtime_job(mock_wildcards, attempt=2)
        medium_result = self.medium_runtime_job(mock_wildcards, attempt=2)
        high_result = self.high_runtime_job(mock_wildcards, attempt=2)

        assert low_result == 4  # 2 * 2
        assert medium_result == 20  # 2 * 10
        assert high_result == 60  # 2 * 30

        # Test that result can be used as runtime allocation (positive integer)
        assert low_result > 0
        assert medium_result > 0
        assert high_result > 0

        assert isinstance(low_result, int)
        assert isinstance(medium_result, int)
        assert isinstance(high_result, int)
