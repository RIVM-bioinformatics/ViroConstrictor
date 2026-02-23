import configparser
import json
import os
import subprocess
import sys
from unittest.mock import MagicMock, Mock, mock_open, patch
from urllib.error import URLError

import packaging.version
import pytest

from ViroConstrictor.update import fetch_online_metadata, post_install, update


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
        mock_urlopen.side_effect = URLError("Connection failed")

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

    @patch("ViroConstrictor.update.subprocess.run")
    @patch("ViroConstrictor.update.sys.exit")
    @patch("ViroConstrictor.update.print")
    def test_post_install_execution(self, mock_print, mock_exit, mock_subprocess):
        """Test post-install execution flow."""
        sysargs = ["viroconstrictor", "--input", "data/"]
        version = packaging.version.parse("2.0.0")

        post_install(sysargs, version)

        mock_print.assert_called_once()
        assert "ViroConstrictor updated to version" in mock_print.call_args[0][0]
        assert "2.0.0" in mock_print.call_args[0][0]
        mock_subprocess.assert_called_once_with(sysargs)
        mock_exit.assert_called_once_with(0)


class TestUpdate:
    """Test the main update function."""

    def create_mock_config(self, auto_update="no", ask_for_update="yes"):
        """Create a mock configuration object."""
        config = Mock(spec=configparser.ConfigParser)
        config.__getitem__ = Mock(return_value={"auto_update": auto_update, "ask_for_update": ask_for_update})
        return config

    def create_mock_metadata(self, version="2.0.0"):
        """Create mock online metadata."""
        return {"distributions": [{"version": version}]}

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_update_auto_continue_no_update_needed(self, mock_fetch):
        """Test auto-update when no update is needed."""
        config = self.create_mock_config(auto_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("1.0.0")

        update(["viroconstrictor"], config)

        mock_fetch.assert_called_once()

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    @patch("ViroConstrictor.update.mamba_install")
    @patch("ViroConstrictor.update.post_install")
    @patch("ViroConstrictor.update.log")
    def test_update_auto_continue_successful_update(self, mock_log, mock_post_install, mock_mamba, mock_fetch):
        """Test successful auto-update."""
        config = self.create_mock_config(auto_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")
        mock_mamba.return_value = True

        sysargs = ["viroconstrictor", "--input", "data/"]
        update(sysargs, config)

        mock_fetch.assert_called_once()
        mock_log.info.assert_called()
        mock_mamba.assert_called_once()
        mock_post_install.assert_called_once_with(sysargs, packaging.version.parse("2.0.0"))

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    @patch("ViroConstrictor.update.mamba_install")
    @patch("ViroConstrictor.update.log")
    def test_update_auto_continue_failed_update(self, mock_log, mock_mamba, mock_fetch):
        """Test failed auto-update."""
        config = self.create_mock_config(auto_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")
        mock_mamba.return_value = False

        update(["viroconstrictor"], config)

        mock_fetch.assert_called_once()
        mock_mamba.assert_called_once()
        # Should not call post_install on failure

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    @patch("ViroConstrictor.update.mamba_install")
    @patch("ViroConstrictor.update.log")
    def test_update_auto_continue_mamba_exception(self, mock_log, mock_mamba, mock_fetch):
        """Test auto-update with mamba exception."""
        config = self.create_mock_config(auto_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")
        mock_mamba.side_effect = Exception("Mamba failed")

        update(["viroconstrictor"], config)

        mock_fetch.assert_called_once()
        mock_mamba.assert_called_once()
        mock_log.error.assert_called()
        mock_log.warning.assert_called()

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_update_auto_continue_no_metadata(self, mock_fetch):
        """Test auto-update when metadata fetch fails."""
        config = self.create_mock_config(auto_update="yes")
        mock_fetch.return_value = None

        result = update(["viroconstrictor"], config)

        assert result is None
        mock_fetch.assert_called_once()

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    def test_update_no_prompt_early_return(self):
        """Test early return when auto_update=False and ask_for_update=False."""
        config = self.create_mock_config(auto_update="no", ask_for_update="no")

        result = update(["viroconstrictor"], config)

        assert result is None

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    @patch("ViroConstrictor.update.AskPrompts")
    @patch("ViroConstrictor.update.mamba_install")
    @patch("ViroConstrictor.update.post_install")
    @patch("ViroConstrictor.update.log")
    def test_update_prompt_user_accepts(self, mock_log, mock_post_install, mock_mamba, mock_ask, mock_fetch):
        """Test prompted update when user accepts."""
        config = self.create_mock_config(auto_update="no", ask_for_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")
        mock_ask.return_value = "yes"
        mock_mamba.return_value = True

        sysargs = ["viroconstrictor", "--input", "data/"]
        update(sysargs, config)

        mock_fetch.assert_called_once()
        mock_ask.assert_called_once()
        mock_log.info.assert_called()
        mock_mamba.assert_called_once()
        mock_post_install.assert_called_once_with(sysargs, packaging.version.parse("2.0.0"))

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    @patch("ViroConstrictor.update.AskPrompts")
    @patch("ViroConstrictor.update.log")
    def test_update_prompt_user_declines(self, mock_log, mock_ask, mock_fetch):
        """Test prompted update when user declines."""
        config = self.create_mock_config(auto_update="no", ask_for_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")
        mock_ask.return_value = "no"

        update(["viroconstrictor"], config)

        mock_fetch.assert_called_once()
        mock_ask.assert_called_once()
        mock_log.info.assert_called_with("Skipping update to version: [bold yellow]2.0.0[/bold yellow]")

    @patch("ViroConstrictor.update.__version__", "2.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_update_prompt_no_update_needed(self, mock_fetch):
        """Test prompted update when no update is needed."""
        config = self.create_mock_config(auto_update="no", ask_for_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")

        result = update(["viroconstrictor"], config)

        assert result is None
        mock_fetch.assert_called_once()

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_update_prompt_no_metadata(self, mock_fetch):
        """Test prompted update when metadata fetch fails."""
        config = self.create_mock_config(auto_update="no", ask_for_update="yes")
        mock_fetch.return_value = None

        result = update(["viroconstrictor"], config)

        assert result is None
        mock_fetch.assert_called_once()

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    @patch("ViroConstrictor.update.AskPrompts")
    @patch("ViroConstrictor.update.mamba_install")
    @patch("ViroConstrictor.update.log")
    def test_update_prompt_mamba_exception(self, mock_log, mock_mamba, mock_ask, mock_fetch):
        """Test prompted update with mamba exception."""
        config = self.create_mock_config(auto_update="no", ask_for_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")
        mock_ask.return_value = "yes"
        mock_mamba.side_effect = Exception("Mamba failed")

        update(["viroconstrictor"], config)

        mock_fetch.assert_called_once()
        mock_ask.assert_called_once()
        mock_mamba.assert_called_once()
        mock_log.error.assert_called()
        mock_log.warning.assert_called()

    def test_version_comparison_logic(self):
        """Test that packaging.version comparison works correctly."""
        # Test basic version comparison
        v1 = packaging.version.parse("1.0.0")
        v2 = packaging.version.parse("2.0.0")
        v3 = packaging.version.parse("1.0.1")

        assert v1 < v2
        assert v1 < v3
        assert v2 > v1
        assert v3 > v1
        assert v1 == packaging.version.parse("1.0.0")

    @patch.dict(os.environ, {"CONDA_PREFIX": "/opt/conda"})
    @patch("ViroConstrictor.update.__prog__", "ViroConstrictor")
    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.fetch_online_metadata")
    @patch("ViroConstrictor.update.mamba_install")
    @patch("ViroConstrictor.update.post_install")
    def test_mamba_install_parameters(self, mock_post_install, mock_mamba, mock_fetch):
        """Test that mamba_install is called with correct parameters."""
        config = self.create_mock_config(auto_update="yes")
        mock_fetch.return_value = self.create_mock_metadata("2.0.0")
        mock_mamba.return_value = True

        update(["viroconstrictor"], config)

        mock_mamba.assert_called_once_with(
            "/opt/conda",
            ("viroconstrictor 2.0.0",),
            ("bioconda", "conda-forge"),
        )


class TestUpdateIntegration:
    """Integration tests for the update module."""

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.request.urlopen")
    @patch("ViroConstrictor.update.mamba_install")
    @patch("ViroConstrictor.update.post_install")
    def test_full_auto_update_flow(self, mock_post_install, mock_mamba, mock_urlopen):
        """Test complete auto-update flow with real-like data."""
        # Mock successful API response
        mock_response = Mock()
        mock_response.read.return_value = json.dumps({"distributions": [{"version": "2.1.0"}]}).encode()
        mock_urlopen.return_value = mock_response
        mock_mamba.return_value = True

        config = Mock(spec=configparser.ConfigParser)
        config.__getitem__ = Mock(return_value={"auto_update": "yes", "ask_for_update": "yes"})

        sysargs = ["viroconstrictor", "--input", "test_data/"]
        update(sysargs, config)

        # Verify the complete flow
        mock_urlopen.assert_called_once()
        mock_mamba.assert_called_once()
        mock_post_install.assert_called_once_with(sysargs, packaging.version.parse("2.1.0"))

    def test_config_parsing_variations(self):
        """Test different configuration parsing scenarios."""
        # Test case sensitivity
        config1 = Mock(spec=configparser.ConfigParser)
        config1.__getitem__ = Mock(return_value={"auto_update": "YES", "ask_for_update": "NO"})

        config2 = Mock(spec=configparser.ConfigParser)
        config2.__getitem__ = Mock(return_value={"auto_update": "false", "ask_for_update": "true"})

        # These should not raise exceptions when called
        with patch("ViroConstrictor.update.fetch_online_metadata", return_value=None):
            update(["test"], config1)
            update(["test"], config2)
