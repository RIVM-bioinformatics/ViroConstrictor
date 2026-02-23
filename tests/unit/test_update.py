import configparser
import json
import os
from unittest.mock import Mock, patch

import packaging.version
import pytest

from ViroConstrictor.update import (
    _build_update_cmd,
    _get_online_version,
    _handle_update_result,
    _prompt_for_update,
    _run_update,
    _update_if_available,
    fetch_online_metadata,
    post_install,
    update,
)


def _build_config(auto_update: str = "no", ask_for_update: str = "no") -> configparser.ConfigParser:
    conf = configparser.ConfigParser()
    conf["GENERAL"] = {
        "auto_update": auto_update,
        "ask_for_update": ask_for_update,
    }
    return conf


class TestFetchOnlineMetadata:
    """Test the fetch_online_metadata function."""

    @patch("ViroConstrictor.update.request.urlopen")
    def test_fetch_online_metadata_success(self, mock_urlopen):
        """Test successful metadata fetch."""
        mock_response = Mock()
        mock_response.read.return_value = b'{"version": "1.0.0", "distributions": []}'
        mock_urlopen.return_value = mock_response

        result = fetch_online_metadata()

        assert result == {"version": "1.0.0", "distributions": []}
        mock_urlopen.assert_called_once()

    @patch("ViroConstrictor.update.request.urlopen")
    @patch("ViroConstrictor.update.log")
    def test_fetch_online_metadata_connection_error(self, mock_log, mock_urlopen):
        """Test metadata fetch with connection error."""
        mock_urlopen.side_effect = Exception("Connection failed")

        result = fetch_online_metadata()

        assert result is None
        mock_log.warning.assert_called_once()
        assert "Unable to connect to Anaconda API" in mock_log.warning.call_args[0][0]

    @patch("ViroConstrictor.update.request.urlopen")
    @patch("ViroConstrictor.update.log")
    def test_fetch_online_metadata_json_error(self, mock_log, mock_urlopen):
        """Test metadata fetch with invalid JSON."""
        mock_response = Mock()
        mock_response.read.return_value = b"invalid json"
        mock_urlopen.return_value = mock_response

        with pytest.raises(json.JSONDecodeError):
            result = fetch_online_metadata()


class TestPostInstall:
    """Test the post_install function."""

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update.os.execv")
    @patch("ViroConstrictor.update.shutil.which")
    @patch("ViroConstrictor.update.sys.exit")
    def test_post_install_execution(self, mock_exit, mock_which, mock_execv, mock_log):
        """Test post-install execution flow."""
        sysargs = ["viroconstrictor", "--input", "data/"]
        version = packaging.version.parse("2.0.0")
        mock_which.return_value = "/opt/conda/bin/viroconstrictor"

        post_install(sysargs, version)

        mock_log.info.assert_called_once_with("ViroConstrictor updated to version [bold yellow]2.0.0[/bold yellow]")
        mock_execv.assert_called_once_with(
            "/opt/conda/bin/viroconstrictor",
            ["/opt/conda/bin/viroconstrictor", "--input", "data/"],
        )
        mock_exit.assert_not_called()

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update.os.execv")
    @patch("ViroConstrictor.update.shutil.which")
    def test_post_install_missing_binary(self, mock_which, mock_execv, mock_log):
        """Test post-install exits when executable is missing."""
        mock_which.return_value = None

        with pytest.raises(SystemExit):
            post_install(["viroconstrictor"], packaging.version.parse("2.0.0"))

        mock_log.error.assert_called_once_with("Could not find viroconstrictor executable after update")
        mock_execv.assert_not_called()


class TestHelperFunctions:
    """Test update helper functions."""

    @patch("ViroConstrictor.update.shutil.which")
    def test_build_update_cmd_prefers_mamba(self, mock_which):
        """Test update command prefers mamba when available."""
        mock_which.return_value = "/usr/bin/mamba"

        cmd = _build_update_cmd(packaging.version.parse("2.0.0"), "/opt/conda")

        assert cmd[0] == "mamba"
        assert "viroconstrictor=2.0.0" in cmd

    @patch("ViroConstrictor.update.shutil.which")
    def test_build_update_cmd_falls_back_to_conda(self, mock_which):
        """Test update command falls back to conda when mamba is unavailable."""
        mock_which.return_value = None

        cmd = _build_update_cmd(packaging.version.parse("2.0.0"), "/opt/conda")

        assert cmd[0] == "conda"

    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_get_online_version_returns_version(self, mock_fetch):
        """Test extracting online version from metadata."""
        mock_fetch.return_value = {"distributions": [{"version": "2.1.0"}]}

        online_version = _get_online_version()

        assert online_version == packaging.version.parse("2.1.0")

    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_get_online_version_no_metadata(self, mock_fetch):
        """Test online version extraction when metadata is unavailable."""
        mock_fetch.return_value = None

        online_version = _get_online_version()

        assert online_version is None

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update.subprocess.run")
    @patch.dict(os.environ, {}, clear=True)
    def test_run_update_requires_conda_prefix(self, mock_subprocess, mock_log):
        """Test update run fails gracefully without CONDA_PREFIX."""
        result = _run_update(packaging.version.parse("2.0.0"))

        assert result is False
        assert mock_subprocess.call_count == 1  # pip uninstall only
        mock_log.error.assert_called_once()

    @patch("ViroConstrictor.update.subprocess.run")
    @patch("ViroConstrictor.update._build_update_cmd")
    @patch.dict(os.environ, {"CONDA_PREFIX": "/opt/conda"}, clear=True)
    def test_run_update_success(self, mock_build_cmd, mock_subprocess):
        """Test successful update command execution."""
        mock_build_cmd.return_value = ["mamba", "install"]
        mock_subprocess.side_effect = [
            Mock(returncode=0),
            Mock(returncode=0, stderr=""),
        ]

        result = _run_update(packaging.version.parse("2.0.0"))

        assert result is True
        assert mock_subprocess.call_count == 2

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update.subprocess.run")
    @patch("ViroConstrictor.update._build_update_cmd")
    @patch.dict(os.environ, {"CONDA_PREFIX": "/opt/conda"}, clear=True)
    def test_run_update_failure_logs_stderr(self, mock_build_cmd, mock_subprocess, mock_log):
        """Test failed update command logs stderr."""
        mock_build_cmd.return_value = ["mamba", "install"]
        mock_subprocess.side_effect = [
            Mock(returncode=0),
            Mock(returncode=1, stderr="solver failed"),
        ]

        result = _run_update(packaging.version.parse("2.0.0"))

        assert result is False
        mock_log.error.assert_called_once()

    @patch("ViroConstrictor.update._get_online_version")
    @patch("ViroConstrictor.update.AskPrompts")
    def test_prompt_for_update_accepts(self, mock_ask, mock_online_version):
        """Test prompting user for update acceptance."""
        mock_online_version.return_value = packaging.version.parse("2.0.0")
        mock_ask.return_value = "yes"

        result = _prompt_for_update(packaging.version.parse("1.0.0"))

        assert result == packaging.version.parse("2.0.0")

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update._get_online_version")
    @patch("ViroConstrictor.update.AskPrompts")
    def test_prompt_for_update_declines(self, mock_ask, mock_online_version, mock_log):
        """Test prompting user for update decline."""
        mock_online_version.return_value = packaging.version.parse("2.0.0")
        mock_ask.return_value = "no"

        result = _prompt_for_update(packaging.version.parse("1.0.0"))

        assert result is None
        mock_log.info.assert_called_once_with("Skipping update to version: [bold yellow]2.0.0[/bold yellow]")


class TestUpdate:
    """Test the main update function."""

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update._update_if_available")
    def test_update_auto_continue_calls_update_if_available(self, mock_update_if_available):
        """Test auto-update delegates to _update_if_available."""
        config = _build_config(auto_update="yes", ask_for_update="yes")

        update(["viroconstrictor"], config)

        mock_update_if_available.assert_called_once_with(["viroconstrictor"], packaging.version.parse("1.0.0"))

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    def test_update_no_prompt_early_return(self):
        """Test early return when auto_update=False and ask_for_update=False."""
        config = _build_config(auto_update="no", ask_for_update="no")

        result = update(["viroconstrictor"], config)

        assert result is None

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update._prompt_for_update")
    @patch("ViroConstrictor.update._run_update")
    @patch("ViroConstrictor.update._handle_update_result")
    @patch("ViroConstrictor.update.log")
    def test_update_prompt_user_accepts(self, mock_log, mock_handle_update_result, mock_run_update, mock_prompt_for_update):
        """Test prompted update when user accepts."""
        config = _build_config(auto_update="no", ask_for_update="yes")
        mock_prompt_for_update.return_value = packaging.version.parse("2.0.0")
        mock_run_update.return_value = True

        sysargs = ["viroconstrictor", "--input", "data/"]
        update(sysargs, config)

        mock_prompt_for_update.assert_called_once_with(packaging.version.parse("1.0.0"))
        mock_run_update.assert_called_once_with(packaging.version.parse("2.0.0"))
        mock_handle_update_result.assert_called_once_with(sysargs, packaging.version.parse("2.0.0"), True)
        mock_log.info.assert_called_once_with("Updating ViroConstrictor to latest version: [bold yellow]2.0.0[/bold yellow]")

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update._prompt_for_update")
    @patch("ViroConstrictor.update._run_update")
    def test_update_prompt_user_declines(self, mock_run_update, mock_prompt_for_update):
        """Test prompted update when user declines."""
        config = _build_config(auto_update="no", ask_for_update="yes")
        mock_prompt_for_update.return_value = None

        update(["viroconstrictor"], config)

        mock_prompt_for_update.assert_called_once()
        mock_run_update.assert_not_called()

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update._prompt_for_update")
    @patch("ViroConstrictor.update._run_update")
    @patch("ViroConstrictor.update.log")
    def test_update_prompt_run_exception(self, mock_log, mock_run_update, mock_prompt_for_update):
        """Test prompted update with solver exception."""
        config = _build_config(auto_update="no", ask_for_update="yes")
        mock_prompt_for_update.return_value = packaging.version.parse("2.0.0")
        mock_run_update.side_effect = Exception("solver failed")

        update(["viroconstrictor"], config)

        mock_prompt_for_update.assert_called_once()
        mock_run_update.assert_called_once_with(packaging.version.parse("2.0.0"))
        mock_log.error.assert_called()
        mock_log.warning.assert_called()

    @patch("ViroConstrictor.update.post_install")
    @patch("ViroConstrictor.update.log")
    def test_handle_update_result_failure_logs_error(self, mock_log, mock_post_install):
        """Test handling failed update result."""
        _handle_update_result(["viroconstrictor"], packaging.version.parse("2.0.0"), False)

        mock_post_install.assert_not_called()
        mock_log.error.assert_called_once()


class TestUpdateIntegration:
    """Integration tests for the update module."""

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.request.urlopen")
    @patch("ViroConstrictor.update._run_update")
    @patch("ViroConstrictor.update.post_install")
    def test_full_auto_update_flow(self, mock_post_install, mock_run_update, mock_urlopen):
        """Test complete auto-update flow with real-like data."""
        # Mock successful API response
        mock_response = Mock()
        mock_response.read.return_value = json.dumps({"distributions": [{"version": "2.1.0"}]}).encode()
        mock_urlopen.return_value = mock_response
        mock_run_update.return_value = True

        config = _build_config(auto_update="yes", ask_for_update="yes")

        sysargs = ["viroconstrictor", "--input", "test_data/"]
        update(sysargs, config)

        # Verify the complete flow
        mock_urlopen.assert_called_once()
        mock_run_update.assert_called_once_with(packaging.version.parse("2.1.0"))
        mock_post_install.assert_called_once_with(sysargs, packaging.version.parse("2.1.0"))

    def test_config_parsing_variations(self):
        """Test different configuration parsing scenarios."""
        # Test case and value variations accepted by ConfigParser.getboolean
        config1 = _build_config(auto_update="YES", ask_for_update="NO")
        config2 = _build_config(auto_update="false", ask_for_update="true")

        # These should not raise exceptions when called
        with patch("ViroConstrictor.update._update_if_available"):
            update(["test"], config1)

        with patch("ViroConstrictor.update._prompt_for_update", return_value=None):
            update(["test"], config2)


class TestUpdateIfAvailable:
    """Tests for the internal _update_if_available workflow."""

    @patch("ViroConstrictor.update._get_online_version")
    @patch("ViroConstrictor.update._run_update")
    @patch("ViroConstrictor.update.post_install")
    def test_update_if_available_success(self, mock_post_install, mock_run_update, mock_get_online):
        """Test successful update path when online version is newer."""
        mock_get_online.return_value = packaging.version.parse("2.0.0")
        mock_run_update.return_value = True

        _update_if_available(["viroconstrictor", "--help"], packaging.version.parse("1.0.0"))

        mock_run_update.assert_called_once_with(packaging.version.parse("2.0.0"))
        mock_post_install.assert_called_once_with(["viroconstrictor", "--help"], packaging.version.parse("2.0.0"))

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update._get_online_version")
    @patch("ViroConstrictor.update._run_update")
    def test_update_if_available_exception(self, mock_run_update, mock_get_online, mock_log):
        """Test exception handling in _update_if_available."""
        mock_get_online.return_value = packaging.version.parse("2.0.0")
        mock_run_update.side_effect = Exception("solver failed")

        _update_if_available(["viroconstrictor"], packaging.version.parse("1.0.0"))

        mock_log.error.assert_called_once()
        mock_log.warning.assert_called_once()
