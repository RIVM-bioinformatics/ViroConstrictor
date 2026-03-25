"""Unit tests for run report PDF helpers and report generation."""

import configparser
import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest
from fpdf import FPDF
from snakemake.api import OutputSettings, ResourceSettings

from ViroConstrictor.parser import CLIparser
from ViroConstrictor.runreport import PDF, WriteReport, analysis_details, directory_sections


class TestPDFClass:
    """Validate the custom PDF wrapper used for run report generation."""

    def test_pdf_initialization(self):
        """Verify that the PDF wrapper initializes as an FPDF instance."""
        pdf = PDF()
        assert isinstance(pdf, FPDF)
        assert hasattr(pdf, "timestamp")
        assert hasattr(pdf, "header")

    def test_timestamp_format(self):
        """Verify that the timestamp uses the expected report format."""
        pdf = PDF()
        timestamp = pdf.timestamp()

        # Should match format: "DD-MM-YYYY HH:MM"
        assert len(timestamp) == 16
        assert timestamp[2] == "-"
        assert timestamp[5] == "-"
        assert timestamp[10] == " "
        assert timestamp[13] == ":"

        # Verify it's a valid datetime by parsing it back
        datetime.strptime(timestamp, "%d-%m-%Y %H:%M")

    @patch("ViroConstrictor.runreport.__version__", "2.0.0")
    def test_header_content(self):
        """Verify that the header renders the report title and version."""
        pdf = PDF()

        # Mock the PDF methods to capture calls
        pdf.set_font = Mock()
        pdf.cell = Mock()

        pdf.header()

        # Verify set_font was called multiple times with different styles
        assert pdf.set_font.call_count >= 4

        # Verify cell calls include expected content
        cell_calls = [call[0] for call in pdf.cell.call_args_list]

        # Check for ViroConstrictor title
        title_calls = [call for call in cell_calls if len(call) > 2 and "ViroConstrictor" in str(call[2])]
        assert len(title_calls) > 0

        # Check for version info
        version_calls = [call for call in cell_calls if len(call) > 2 and "2.0.0" in str(call[2])]
        assert len(version_calls) > 0


class TestDirectorySections:
    """Validate directory section formatting in the report body."""

    def create_mock_pdf(self):
        """Return a PDF mock with the methods used by directory sections."""
        pdf = Mock(spec=PDF)
        pdf.set_font = Mock()
        pdf.cell = Mock()
        return pdf

    def test_directory_sections_complete_data(self):
        """Verify that directory sections render when data is complete."""
        pdf = self.create_mock_pdf()
        contents = {0: ["Output directory:", "Description", "/path/to/output"], 1: ["Input directory:", "Description", "/path/to/input"]}

        result_pdf = directory_sections(pdf, 0, contents)

        assert result_pdf == pdf
        assert pdf.set_font.called
        assert pdf.cell.called

    def test_directory_sections_missing_data(self):
        """Verify that missing section data is ignored safely."""
        pdf = self.create_mock_pdf()
        contents = {1: ["Input directory:", "Description", "/path/to/input"]}

        # Should not crash when iteration 0 is missing
        result_pdf = directory_sections(pdf, 0, contents)
        assert result_pdf == pdf

    def test_directory_sections_empty_contents(self):
        """Verify that empty section content leaves the PDF unchanged."""
        pdf = self.create_mock_pdf()
        contents = {}

        result_pdf = directory_sections(pdf, 0, contents)
        assert result_pdf == pdf


class TestAnalysisDetails:
    """Validate analysis detail formatting in the report body."""

    def create_mock_pdf(self):
        """Return a PDF mock with the methods used by analysis details."""
        pdf = Mock(spec=PDF)
        pdf.set_font = Mock()
        pdf.cell = Mock()
        return pdf

    def test_analysis_details_basic(self):
        """Verify that analysis details render header and body text."""
        pdf = self.create_mock_pdf()
        header = "Test Header:"
        text = "Test content"

        result_pdf = analysis_details(pdf, header, text)

        assert result_pdf == pdf
        assert pdf.set_font.call_count == 2  # Bold for header, normal for text
        assert pdf.cell.call_count == 2  # One for header, one for text

    def test_analysis_details_empty_text(self):
        """Verify that empty analysis text is handled without error."""
        pdf = self.create_mock_pdf()

        result_pdf = analysis_details(pdf, "Header:", "")
        assert result_pdf == pdf

    def test_analysis_details_special_characters(self):
        """Verify that path-like text is rendered without modification."""
        pdf = self.create_mock_pdf()

        result_pdf = analysis_details(pdf, "Path:", "/path/with/special-chars_123")
        assert result_pdf == pdf


class TestWriteReport:
    """Validate the top-level report writer orchestration."""

    def create_mock_config(self):
        """Return a configuration mock with the expected sections."""
        config = Mock(spec=configparser.ConfigParser)
        config.__getitem__ = Mock(side_effect=lambda key: {"COMPUTING": {"compmode": "local", "queuename": "default"}}[key])
        return config

    def create_mock_resource_settings(self):
        """Return a Snakemake resource settings mock."""
        settings = Mock(spec=ResourceSettings)
        settings.cores = 4
        return settings

    def create_mock_output_settings(self, dryrun=False):
        """Return a Snakemake output settings mock."""
        settings = Mock(spec=OutputSettings)
        settings.dryrun = dryrun
        return settings

    def create_mock_cli_parser(self):
        """Return a CLI parser mock with the flags used by the report."""
        parser = Mock(spec=CLIparser)
        parser.flags = Mock()
        parser.flags.platform = "illumina"
        parser.flags.amplicon_type = "V4"
        return parser

    @patch("ViroConstrictor.runreport.os.chdir")
    @patch("ViroConstrictor.runreport.os.getcwd")
    @patch("ViroConstrictor.runreport.PDF")
    @patch("ViroConstrictor.runreport.sys.argv", ["viroconstrictor", "--input", "data/", "--output", "results/"])
    def test_write_report_success_local(self, mock_pdf_class, mock_getcwd, mock_chdir):
        """Verify that local report generation writes the expected PDF."""
        # Setup mocks
        mock_getcwd.return_value = "/current/dir"
        mock_pdf = Mock()
        mock_pdf_class.return_value = mock_pdf

        workingdir = "/output/dir"
        inpath = "/input/dir"
        startpath = "/start/dir"
        config = self.create_mock_config()
        resource_settings = self.create_mock_resource_settings()
        output_settings = self.create_mock_output_settings()
        cli_parser = self.create_mock_cli_parser()

        # Execute
        WriteReport(
            workingdir=workingdir,
            inpath=inpath,
            startpath=startpath,
            conf=config,
            snakemake_resource_settings=resource_settings,
            snakemake_output_conf=output_settings,
            inputs_config=cli_parser,
            status="Success",
        )

        # Verify directory change
        mock_chdir.assert_called_once_with(workingdir)

        # Verify PDF creation and setup
        mock_pdf_class.assert_called_once()
        mock_pdf.set_author.assert_called_once_with("SARS2seq")
        mock_pdf.add_page.assert_called_once()
        mock_pdf.output.assert_called_once_with(name="Runinfomation.pdf")

    @patch("ViroConstrictor.runreport.os.chdir")
    @patch("ViroConstrictor.runreport.os.getcwd")
    @patch("ViroConstrictor.runreport.PDF")
    def test_write_report_dryrun(self, mock_pdf_class, mock_getcwd, mock_chdir):
        """Verify that dry-run execution skips the directory change."""
        mock_getcwd.return_value = "/output/dir"
        mock_pdf = Mock()
        mock_pdf_class.return_value = mock_pdf

        config = self.create_mock_config()
        resource_settings = self.create_mock_resource_settings()
        output_settings = self.create_mock_output_settings(dryrun=True)
        cli_parser = self.create_mock_cli_parser()

        WriteReport(
            workingdir="/output/dir",
            inpath="/input/dir",
            startpath="/start/dir",
            conf=config,
            snakemake_resource_settings=resource_settings,
            snakemake_output_conf=output_settings,
            inputs_config=cli_parser,
            status="Success",
        )

        # Should not change directory if already in workingdir
        mock_chdir.assert_not_called()
        mock_pdf.output.assert_called_once_with(name="Runinfomation.pdf")

    @patch("ViroConstrictor.runreport.os.chdir")
    @patch("ViroConstrictor.runreport.os.getcwd")
    @patch("ViroConstrictor.runreport.PDF")
    def test_write_report_grid_execution(self, mock_pdf_class, mock_getcwd, mock_chdir):
        """Verify that grid execution metadata is accepted."""
        mock_getcwd.return_value = "/current/dir"
        mock_pdf = Mock()
        mock_pdf_class.return_value = mock_pdf

        # Grid execution config
        config = Mock(spec=configparser.ConfigParser)
        config.__getitem__ = Mock(side_effect=lambda key: {"COMPUTING": {"compmode": "grid", "queuename": "gpu_queue"}}[key])

        resource_settings = self.create_mock_resource_settings()
        output_settings = self.create_mock_output_settings()
        cli_parser = self.create_mock_cli_parser()

        WriteReport(
            workingdir="/output/dir",
            inpath="/input/dir",
            startpath="/start/dir",
            conf=config,
            snakemake_resource_settings=resource_settings,
            snakemake_output_conf=output_settings,
            inputs_config=cli_parser,
            status="Success",
        )

        mock_pdf.output.assert_called_once_with(name="Runinfomation.pdf")

    @patch("ViroConstrictor.runreport.os.chdir")
    @patch("ViroConstrictor.runreport.os.getcwd")
    @patch("ViroConstrictor.runreport.PDF")
    def test_write_report_failed_status(self, mock_pdf_class, mock_getcwd, mock_chdir):
        """Verify that failed status still produces a report."""
        mock_getcwd.return_value = "/output/dir"
        mock_pdf = Mock()
        mock_pdf_class.return_value = mock_pdf

        config = self.create_mock_config()
        resource_settings = self.create_mock_resource_settings()
        output_settings = self.create_mock_output_settings()
        cli_parser = self.create_mock_cli_parser()

        WriteReport(
            workingdir="/output/dir",
            inpath="/input/dir",
            startpath="/start/dir",
            conf=config,
            snakemake_resource_settings=resource_settings,
            snakemake_output_conf=output_settings,
            inputs_config=cli_parser,
            status="Failed",
        )

        mock_pdf.output.assert_called_once_with(name="Runinfomation.pdf")

    @patch("ViroConstrictor.runreport.os.chdir")
    @patch("ViroConstrictor.runreport.os.getcwd")
    @patch("ViroConstrictor.runreport.PDF")
    @patch("ViroConstrictor.runreport.sys.argv", ["viroconstrictor", "--input", "test data/", "--output", "results/", "--amplicon-type", "v3"])
    def test_write_report_command_capture(self, mock_pdf_class, mock_getcwd, mock_chdir):
        """Verify that the command line is captured in the report."""
        mock_getcwd.return_value = "/output/dir"
        mock_pdf = Mock()
        mock_pdf_class.return_value = mock_pdf

        config = self.create_mock_config()
        resource_settings = self.create_mock_resource_settings()
        output_settings = self.create_mock_output_settings()
        cli_parser = self.create_mock_cli_parser()

        WriteReport(
            workingdir="/output/dir",
            inpath="/input/dir",
            startpath="/start/dir",
            conf=config,
            snakemake_resource_settings=resource_settings,
            snakemake_output_conf=output_settings,
            inputs_config=cli_parser,
            status="Success",
        )

        # Verify multi_cell was called (for command display)
        mock_pdf.multi_cell.assert_called_once()

        # Get the command that was written
        command_call = mock_pdf.multi_cell.call_args[0][2]
        assert "viroconstrictor" in command_call
        assert "--input" in command_call
        assert "test data/" in command_call


class TestReportIntegration:
    """Validate end-to-end report generation with realistic inputs."""

    def test_report_generation_integration(self):
        """Verify that a report file is created for realistic inputs."""
        with tempfile.TemporaryDirectory() as temp_dir:
            workdir = Path(temp_dir) / "output"
            workdir.mkdir()

            # Create realistic test data
            config = configparser.ConfigParser()
            config.add_section("COMPUTING")
            config.set("COMPUTING", "compmode", "local")
            config.set("COMPUTING", "queuename", "default")

            resource_settings = Mock(spec=ResourceSettings)
            resource_settings.cores = 8

            output_settings = Mock(spec=OutputSettings)
            output_settings.dryrun = False

            cli_parser = Mock(spec=CLIparser)
            cli_parser.flags = Mock()
            cli_parser.flags.platform = "nanopore"
            cli_parser.flags.amplicon_type = "V4"

            with patch("ViroConstrictor.runreport.sys.argv", ["viroconstrictor", "--input", "testdata/"]):
                # Should not raise exceptions
                WriteReport(
                    workingdir=str(workdir),
                    inpath="/test/input",
                    startpath="/test/start",
                    conf=config,
                    snakemake_resource_settings=resource_settings,
                    snakemake_output_conf=output_settings,
                    inputs_config=cli_parser,
                    status="Success",
                )

            # Verify PDF was created
            pdf_path = workdir / "Runinfomation.pdf"
            assert pdf_path.exists()
            assert pdf_path.stat().st_size > 0

    @patch("ViroConstrictor.runreport.os.getcwd")
    def test_all_platform_types(self, mock_getcwd):
        """Verify that report generation supports all platform values."""
        platforms = ["illumina", "nanopore", "iontorrent"]

        for platform in platforms:
            with tempfile.TemporaryDirectory() as temp_dir:
                workdir = Path(temp_dir)

                config = configparser.ConfigParser()
                config.add_section("COMPUTING")
                config.set("COMPUTING", "compmode", "local")

                resource_settings = Mock(spec=ResourceSettings)
                resource_settings.cores = 4

                output_settings = Mock(spec=OutputSettings)
                output_settings.dryrun = False

                cli_parser = Mock(spec=CLIparser)
                cli_parser.flags = Mock()
                cli_parser.flags.platform = platform
                cli_parser.flags.amplicon_type = "custom"

                with patch("ViroConstrictor.runreport.sys.argv", ["viroconstrictor"]):
                    WriteReport(
                        workingdir=str(workdir),
                        inpath="/test",
                        startpath="/test",
                        conf=config,
                        snakemake_resource_settings=resource_settings,
                        snakemake_output_conf=output_settings,
                        inputs_config=cli_parser,
                        status="Success",
                    )

                pdf_path = workdir / "Runinfomation.pdf"
                assert pdf_path.exists(), f"PDF not created for platform: {platform}"


class TestReportErrorHandling:
    """Validate error handling during report generation."""

    @patch("ViroConstrictor.runreport.os.getcwd")
    @patch("ViroConstrictor.runreport.os.chdir")
    @patch("ViroConstrictor.runreport.PDF")
    def test_pdf_creation_failure(self, mock_getcwd, mock_pdf_class, mock_chdir):
        """Verify that PDF construction failures are surfaced."""
        mock_pdf_class.side_effect = Exception("PDF creation failed")

        config = configparser.ConfigParser()
        config.add_section("COMPUTING")
        config.set("COMPUTING", "compmode", "local")

        resource_settings = Mock(spec=ResourceSettings)
        resource_settings.cores = 4

        output_settings = Mock(spec=OutputSettings)
        output_settings.dryrun = False

        cli_parser = Mock(spec=CLIparser)
        cli_parser.flags = Mock()
        cli_parser.flags.platform = "illumina"
        cli_parser.flags.amplicon_type = "V4"

        with pytest.raises(Exception, match="PDF creation failed"):
            WriteReport(
                workingdir="/test",
                inpath="/test",
                startpath="/test",
                conf=config,
                snakemake_resource_settings=resource_settings,
                snakemake_output_conf=output_settings,
                inputs_config=cli_parser,
                status="Success",
            )

    @patch("ViroConstrictor.runreport.os.chdir")
    @patch("ViroConstrictor.runreport.os.getcwd")
    def test_directory_change_failure(self, mock_getcwd, mock_chdir):
        """Verify that directory change failures are surfaced."""
        mock_getcwd.return_value = "/current"
        mock_chdir.side_effect = OSError("Permission denied")

        config = configparser.ConfigParser()
        resource_settings = Mock(spec=ResourceSettings)
        output_settings = Mock(spec=OutputSettings)
        cli_parser = Mock(spec=CLIparser)
        cli_parser.flags = Mock()

        with pytest.raises(OSError, match="Permission denied"):
            WriteReport(
                workingdir="/restricted",
                inpath="/test",
                startpath="/test",
                conf=config,
                snakemake_resource_settings=resource_settings,
                snakemake_output_conf=output_settings,
                inputs_config=cli_parser,
                status="Success",
            )
