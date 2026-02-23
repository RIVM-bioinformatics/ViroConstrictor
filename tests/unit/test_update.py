"""
Unit tests for :mod:`ViroConstrictor.update`.

This module contains unit and integration tests for the update helper
functions and higher-level update flow. Tests cover metadata fetching,
post-install behavior, command construction, prompting, and the overall
update orchestration. Network and subprocess interactions are mocked
where appropriate to keep tests hermetic.
"""

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
    """Tests for :func:`ViroConstrictor.update.fetch_online_metadata`.

    Validate that remote metadata is parsed correctly and that
    connection and JSON decoding errors are handled as expected.
    """

    @patch("ViroConstrictor.update.request.urlopen")
    def test_fetch_online_metadata_success(self, mock_urlopen):
        """Return parsed metadata when the API responds successfully.

        Parameters
        ----------
        mock_urlopen : unittest.mock.Mock
            Mock of ``request.urlopen`` that returns a response whose
            ``read()`` method yields a JSON payload.
        """
        mock_response = Mock()
        mock_response.read.return_value = b'{"version": "1.0.0", "distributions": []}'
        mock_urlopen.return_value = mock_response

        result = fetch_online_metadata()

        assert result == {"version": "1.0.0", "distributions": []}
        mock_urlopen.assert_called_once()

    @patch("ViroConstrictor.update.request.urlopen")
    @patch("ViroConstrictor.update.log")
    def test_fetch_online_metadata_connection_error(self, mock_log, mock_urlopen):
        """Return ``None`` and log a warning when connection fails.

        Parameters
        ----------
        mock_log : unittest.mock.Mock
            Mock logger to assert that a warning was emitted.
        mock_urlopen : unittest.mock.Mock
            Mock configured to raise an exception on call.
        """
        mock_urlopen.side_effect = Exception("Connection failed")

        result = fetch_online_metadata()

        assert result is None
        mock_log.warning.assert_called_once()
        assert "Unable to connect to Anaconda API" in mock_log.warning.call_args[0][0]

    @patch("ViroConstrictor.update.request.urlopen")
    @patch("ViroConstrictor.update.log")
    def test_fetch_online_metadata_json_error(self, mock_log, mock_urlopen):
        """Raise ``JSONDecodeError`` when API returns invalid JSON.

        Parameters
        ----------
        mock_log : unittest.mock.Mock
            Mock logger (not asserted in this test).
        mock_urlopen : unittest.mock.Mock
            Mock returning invalid JSON bytes.
        """
        mock_response = Mock()
        mock_response.read.return_value = b"invalid json"
        mock_urlopen.return_value = mock_response

        with pytest.raises(json.JSONDecodeError):
            result = fetch_online_metadata()


class TestPostInstall:
    """Tests for :func:`ViroConstrictor.update.post_install`.

    Verify that the updated executable is executed when present and
    that a missing executable results in an error/exit.
    """

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update.os.execv")
    @patch("ViroConstrictor.update.shutil.which")
    @patch("ViroConstrictor.update.sys.exit")
    def test_post_install_execution(self, mock_exit, mock_which, mock_execv, mock_log):
        """Invoke ``os.execv`` with the updated executable path.

        Parameters
        ----------
        mock_exit : unittest.mock.Mock
            Mock for ``sys.exit`` to ensure it is not called.
        mock_which : unittest.mock.Mock
            Mock for ``shutil.which`` returning the executable path.
        mock_execv : unittest.mock.Mock
            Mock for ``os.execv`` to observe invocation.
        mock_log : unittest.mock.Mock
            Mock logger to assert info logging.
        """
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
        """Raise ``SystemExit`` and log error if executable missing.

        Parameters
        ----------
        mock_which : unittest.mock.Mock
            Mock for ``shutil.which`` returning ``None``.
        mock_execv : unittest.mock.Mock
            Mock for ``os.execv`` (should not be called).
        mock_log : unittest.mock.Mock
            Mock logger to assert error logging.
        """
        mock_which.return_value = None

        with pytest.raises(SystemExit):
            post_install(["viroconstrictor"], packaging.version.parse("2.0.0"))

        mock_log.error.assert_called_once_with("Could not find viroconstrictor executable after update")
        mock_execv.assert_not_called()


class TestHelperFunctions:
    """Tests for helper utilities in the update module.

    This class covers command construction and lightweight helpers used
    by the update runner.
    """

    @patch("ViroConstrictor.update.shutil.which")
    def test_build_update_cmd_prefers_mamba(self, mock_which):
        """Prefer ``mamba`` when available for building the install command.

        Parameters
        ----------
        mock_which : unittest.mock.Mock
            Mock for ``shutil.which`` returning the path to ``mamba``.
        """
        mock_which.return_value = "/usr/bin/mamba"

        cmd = _build_update_cmd(packaging.version.parse("2.0.0"), "/opt/conda")

        assert cmd[0] == "mamba"
        assert "viroconstrictor=2.0.0" in cmd

    @patch("ViroConstrictor.update.shutil.which")
    def test_build_update_cmd_falls_back_to_conda(self, mock_which):
        """Fall back to ``conda`` when ``mamba`` is not present.

        Parameters
        ----------
        mock_which : unittest.mock.Mock
            Mock for ``shutil.which`` returning ``None``.
        """
        mock_which.return_value = None

        cmd = _build_update_cmd(packaging.version.parse("2.0.0"), "/opt/conda")

        assert cmd[0] == "conda"

    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_get_online_version_returns_version(self, mock_fetch):
        """Return a parsed :class:`packaging.version.Version` from metadata.

        Parameters
        ----------
        mock_fetch : unittest.mock.Mock
            Mock of ``fetch_online_metadata`` returning distributions
            including a version string.
        """
        mock_fetch.return_value = {"distributions": [{"version": "2.1.0"}]}

        online_version = _get_online_version()

        assert online_version == packaging.version.parse("2.1.0")

    @patch("ViroConstrictor.update.fetch_online_metadata")
    def test_get_online_version_no_metadata(self, mock_fetch):
        """Return ``None`` when no metadata could be fetched.

        Parameters
        ----------
        mock_fetch : unittest.mock.Mock
            Mock of ``fetch_online_metadata`` returning ``None``.
        """
        mock_fetch.return_value = None

        online_version = _get_online_version()

        assert online_version is None

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update.subprocess.run")
    @patch.dict(os.environ, {}, clear=True)
    def test_run_update_requires_conda_prefix(self, mock_subprocess, mock_log):
        """Return ``False`` and log an error when ``CONDA_PREFIX`` is unset.

        Parameters
        ----------
        mock_subprocess : unittest.mock.Mock
            Mock for ``subprocess.run``.
        mock_log : unittest.mock.Mock
            Mock logger to assert error logging.
        """
        result = _run_update(packaging.version.parse("2.0.0"))

        assert result is False
        assert mock_subprocess.call_count == 1  # pip uninstall only
        mock_log.error.assert_called_once()

    @patch("ViroConstrictor.update.subprocess.run")
    @patch("ViroConstrictor.update._build_update_cmd")
    @patch.dict(os.environ, {"CONDA_PREFIX": "/opt/conda"}, clear=True)
    def test_run_update_success(self, mock_build_cmd, mock_subprocess):
        """Return ``True`` when both pip uninstall and package-manager steps succeed.

        Parameters
        ----------
        mock_build_cmd : unittest.mock.Mock
            Mock for ``_build_update_cmd`` returning package-manager args.
        mock_subprocess : unittest.mock.Mock
            Mock for ``subprocess.run`` returning successful codes.
        """
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
        """Log package-manager stderr and return ``False`` on failure.

        Parameters
        ----------
        mock_build_cmd : unittest.mock.Mock
            Mock for command builder.
        mock_subprocess : unittest.mock.Mock
            Mock for ``subprocess.run`` simulating failure with ``stderr``.
        mock_log : unittest.mock.Mock
            Mock logger to assert an error was emitted.
        """
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
        """Return the online version when the user accepts the prompt.

        Parameters
        ----------
        mock_ask : unittest.mock.Mock
            Mock for prompt helper returning ``"yes"``.
        mock_online_version : unittest.mock.Mock
            Mock for ``_get_online_version`` returning a newer version.
        """
        mock_online_version.return_value = packaging.version.parse("2.0.0")
        mock_ask.return_value = "yes"

        result = _prompt_for_update(packaging.version.parse("1.0.0"))

        assert result == packaging.version.parse("2.0.0")

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update._get_online_version")
    @patch("ViroConstrictor.update.AskPrompts")
    def test_prompt_for_update_declines(self, mock_ask, mock_online_version, mock_log):
        """Return ``None`` and log an informational message when declined.

        Parameters
        ----------
        mock_ask : unittest.mock.Mock
            Mock for prompt helper returning ``"no"``.
        mock_online_version : unittest.mock.Mock
            Mock for ``_get_online_version`` returning a newer version.
        mock_log : unittest.mock.Mock
            Mock logger to assert info logging.
        """
        mock_online_version.return_value = packaging.version.parse("2.0.0")
        mock_ask.return_value = "no"

        result = _prompt_for_update(packaging.version.parse("1.0.0"))

        assert result is None
        mock_log.info.assert_called_once_with("Skipping update to version: [bold yellow]2.0.0[/bold yellow]")


class TestUpdate:
    """Tests for the top-level :func:`ViroConstrictor.update.update`.

    Validate control flow depending on configuration flags and ensure
    helper functions are invoked correctly.
    """

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update._update_if_available")
    def test_update_auto_continue_calls_update_if_available(self, mock_update_if_available):
        """Delegate to ``_update_if_available`` when auto-update enabled.

        Parameters
        ----------
        mock_update_if_available : unittest.mock.Mock
            Mock for the internal update function.
        """
        config = _build_config(auto_update="yes", ask_for_update="yes")

        update(["viroconstrictor"], config)

        mock_update_if_available.assert_called_once_with(["viroconstrictor"], packaging.version.parse("1.0.0"))

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    def test_update_no_prompt_early_return(self):
        """Return ``None`` when both auto-update and prompts are disabled.

        This ensures ``update`` is a no-op for users who opted out.
        """
        config = _build_config(auto_update="no", ask_for_update="no")

        result = update(["viroconstrictor"], config)

        assert result is None

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update._prompt_for_update")
    @patch("ViroConstrictor.update._run_update")
    @patch("ViroConstrictor.update._handle_update_result")
    @patch("ViroConstrictor.update.log")
    def test_update_prompt_user_accepts(self, mock_log, mock_handle_update_result, mock_run_update, mock_prompt_for_update):
        """Run update and handle its result when user accepts prompt.

        Parameters
        ----------
        mock_log : unittest.mock.Mock
            Mock logger used to assert the informational message.
        mock_handle_update_result : unittest.mock.Mock
            Mock for handling the result after update run.
        mock_run_update : unittest.mock.Mock
            Mock that simulates successful update execution.
        mock_prompt_for_update : unittest.mock.Mock
            Mock returning the online version selected for update.
        """
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
        """Do not run the update when the prompt returns ``None``.

        Parameters
        ----------
        mock_run_update : unittest.mock.Mock
            Mock for the update runner (should not be invoked).
        mock_prompt_for_update : unittest.mock.Mock
            Mock returning ``None`` to indicate user declined.
        """
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
        """Log error and warning when the update runner raises an exception.

        Parameters
        ----------
        mock_log : unittest.mock.Mock
            Mock logger to assert error/warning were emitted.
        mock_run_update : unittest.mock.Mock
            Mock that raises to simulate failure.
        mock_prompt_for_update : unittest.mock.Mock
            Mock returning the online version to attempt update.
        """
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
        """Do not call post-install and log an error on failure.

        Parameters
        ----------
        mock_log : unittest.mock.Mock
            Mock logger to assert error logging.
        mock_post_install : unittest.mock.Mock
            Mock for post-install which should not be called on failure.
        """
        _handle_update_result(["viroconstrictor"], packaging.version.parse("2.0.0"), False)

        mock_post_install.assert_not_called()
        mock_log.error.assert_called_once()


class TestUpdateIntegration:
    """Integration-style tests for the update module.

    These tests patch external boundaries but exercise the real
    top-level ``update`` function to validate composed behaviour.
    """

    @patch("ViroConstrictor.update.__version__", "1.0.0")
    @patch("ViroConstrictor.update.request.urlopen")
    @patch("ViroConstrictor.update._run_update")
    @patch("ViroConstrictor.update.post_install")
    def test_full_auto_update_flow(self, mock_post_install, mock_run_update, mock_urlopen):
        """Complete auto-update path using mocked external calls.

        Parameters
        ----------
        mock_post_install : unittest.mock.Mock
            Mock for the post-install handler.
        mock_run_update : unittest.mock.Mock
            Mock simulating a successful update run.
        mock_urlopen : unittest.mock.Mock
            Mock for network call returning metadata with a newer
            distribution.
        """
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
        """Accept multiple string representations for boolean flags.

        Ensure ``update`` handles case variations like "YES"/"NO" and
        truthy/falsey strings without raising.
        """
        # Test case and value variations accepted by ConfigParser.getboolean
        config1 = _build_config(auto_update="YES", ask_for_update="NO")
        config2 = _build_config(auto_update="false", ask_for_update="true")

        # These should not raise exceptions when called
        with patch("ViroConstrictor.update._update_if_available"):
            update(["test"], config1)

        with patch("ViroConstrictor.update._prompt_for_update", return_value=None):
            update(["test"], config2)


class TestUpdateIfAvailable:
    """Tests for :func:`ViroConstrictor.update._update_if_available`.

    Validate the conditional attempt to update when a newer online
    version is detected and verify exception handling during update.
    """

    @patch("ViroConstrictor.update._get_online_version")
    @patch("ViroConstrictor.update._run_update")
    @patch("ViroConstrictor.update.post_install")
    def test_update_if_available_success(self, mock_post_install, mock_run_update, mock_get_online):
        """Run update and post-install when a newer version is available.

        Parameters
        ----------
        mock_post_install : unittest.mock.Mock
            Mock for post-install handler.
        mock_run_update : unittest.mock.Mock
            Mock for the update runner returning success.
        mock_get_online : unittest.mock.Mock
            Mock returning a newer online version.
        """
        mock_get_online.return_value = packaging.version.parse("2.0.0")
        mock_run_update.return_value = True

        _update_if_available(["viroconstrictor", "--help"], packaging.version.parse("1.0.0"))

        mock_run_update.assert_called_once_with(packaging.version.parse("2.0.0"))
        mock_post_install.assert_called_once_with(["viroconstrictor", "--help"], packaging.version.parse("2.0.0"))

    @patch("ViroConstrictor.update.log")
    @patch("ViroConstrictor.update._get_online_version")
    @patch("ViroConstrictor.update._run_update")
    def test_update_if_available_exception(self, mock_run_update, mock_get_online, mock_log):
        """Log error and warning when the update runner raises.

        Parameters
        ----------
        mock_run_update : unittest.mock.Mock
            Mock configured to raise an exception.
        mock_get_online : unittest.mock.Mock
            Mock returning the newer online version.
        mock_log : unittest.mock.Mock
            Mock logger to assert error and warning calls.
        """
        mock_get_online.return_value = packaging.version.parse("2.0.0")
        mock_run_update.side_effect = Exception("solver failed")

        _update_if_available(["viroconstrictor"], packaging.version.parse("1.0.0"))

        mock_log.error.assert_called_once()
        mock_log.warning.assert_called_once()
