"""Unit tests for the viroconstrictor-data compatibility check.

Tests verify that _check_data_compatibility raises SystemExit with an
actionable message when the installed viroconstrictor-data package is not
compatible with the running ViroConstrictor version or is absent entirely.
"""

from unittest.mock import MagicMock, patch

import pytest

import ViroConstrictor.workflow_config as workflow_config


def _make_mock_data_pkg(compatible_range: str) -> MagicMock:
    """Build a minimal viroconstrictor_data mock that returns the given compatible range."""
    mock = MagicMock()
    mock.get_manifest.return_value = {
        "schema_version": "1",
        "compatible_viroconstrictor": compatible_range,
    }
    return mock


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_compatible_version_passes(monkeypatch: pytest.MonkeyPatch) -> None:
    """No exception is raised when the running version falls within the compatible range."""
    monkeypatch.setattr(workflow_config, "__version__", "1.7.0")

    mock_pkg = _make_mock_data_pkg(">=1.7.0,<2.0.0")
    with patch.dict("sys.modules", {"viroconstrictor_data": mock_pkg}):
        # Must not raise
        workflow_config._check_data_compatibility()


# ---------------------------------------------------------------------------
# Incompatible version
# ---------------------------------------------------------------------------


def test_incompatible_version_raises_systemexit(monkeypatch: pytest.MonkeyPatch) -> None:
    """SystemExit is raised when the running version is outside the compatible range."""
    monkeypatch.setattr(workflow_config, "__version__", "2.0.0")

    mock_pkg = _make_mock_data_pkg(">=1.7.0,<2.0.0")
    with patch.dict("sys.modules", {"viroconstrictor_data": mock_pkg}):
        with pytest.raises(SystemExit) as exc_info:
            workflow_config._check_data_compatibility()

    message = str(exc_info.value)
    assert "2.0.0" in message
    assert ">=1.7.0,<2.0.0" in message


# ---------------------------------------------------------------------------
# Missing package
# ---------------------------------------------------------------------------


def test_missing_package_raises_systemexit_with_install_instructions() -> None:
    """SystemExit is raised (not a raw traceback) when viroconstrictor_data is not installed."""
    with patch.dict("sys.modules", {"viroconstrictor_data": None}):
        with pytest.raises(SystemExit) as exc_info:
            workflow_config._check_data_compatibility()

    message = str(exc_info.value)
    assert "viroconstrictor-data is not installed" in message
    assert "pip install" in message


def test_manifest_missing_compatible_key_raises_systemexit() -> None:
    """SystemExit is raised when compatible_viroconstrictor is missing from manifest."""
    mock_pkg = MagicMock()
    mock_pkg.get_manifest.return_value = {"schema_version": "1"}

    with patch.dict("sys.modules", {"viroconstrictor_data": mock_pkg}):
        with pytest.raises(SystemExit) as exc_info:
            workflow_config._check_data_compatibility()

    message = str(exc_info.value)
    assert "compatible_viroconstrictor" in message
    assert "invalid" in message


def test_manifest_invalid_specifier_raises_systemexit() -> None:
    """SystemExit is raised when compatible_viroconstrictor contains an invalid specifier."""
    mock_pkg = _make_mock_data_pkg("not-a-valid-specifier")

    with patch.dict("sys.modules", {"viroconstrictor_data": mock_pkg}):
        with pytest.raises(SystemExit) as exc_info:
            workflow_config._check_data_compatibility()

    message = str(exc_info.value)
    assert "compatible_viroconstrictor" in message
    assert "not valid" in message


def test_manifest_must_be_dictionary() -> None:
    """SystemExit is raised when get_manifest() does not return a dictionary."""
    mock_pkg = MagicMock()
    mock_pkg.get_manifest.return_value = ["not", "a", "dict"]

    with patch.dict("sys.modules", {"viroconstrictor_data": mock_pkg}):
        with pytest.raises(SystemExit) as exc_info:
            workflow_config._check_data_compatibility()

    message = str(exc_info.value)
    assert "manifest must be a dictionary" in message
