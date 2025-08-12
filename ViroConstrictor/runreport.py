import configparser
import os
import sys
from datetime import datetime
from typing import Literal

from fpdf import FPDF
from snakemake.api import OutputSettings, ResourceSettings

from ViroConstrictor import __version__
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.logging import log


class PDF(FPDF):
    log.debug("Python code :: Run report :: Create run report :: creating run information report.")
    def timestamp(self) -> str:
        return datetime.now().strftime("%d-%m-%Y %H:%M")

    def header(self) -> None:
        self.set_font("Helvetica", size=23, style="B")
        self.cell(0, 20, "ViroConstrictor", align="C", ln=1)
        self.set_font("Helvetica", size=12)
        self.cell(0, 5, "Run information summary", align="C", ln=1)
        self.set_font("Helvetica", size=8)
        self.cell(0, 5, self.timestamp(), align="C", ln=1)
        self.set_font("Helvetica", size=10, style="BU")
        self.cell(0, 15, f"Version: {__version__}", align="C", ln=1)


def directory_sections(pdf: PDF, iteration: int, contents: dict[int, list[str]]) -> PDF:
    pdf.set_font("Helvetica", size=12, style="B")
    if outdir := contents.get(iteration):
        pdf.cell(40, 12, outdir[0], align="L")
    pdf.set_font("Helvetica", size=10)
    if indir := contents.get(iteration):
        pdf.cell(0, 12, indir[1], align="L", ln=1)
    pdf.set_font("Helvetica", size=12, style="I")
    if startdir := contents.get(iteration):
        pdf.cell(0, 5, startdir[2], align="L", ln=1)
    # pdf.cell(0, 5, contents.get(iteration)[2], align="L", ln=1)
    return pdf


def analysis_details(pdf: PDF, header: str, text: str) -> PDF:
    pdf.set_font("Helvetica", size=12, style="B")
    pdf.cell(55, 5, header, align="L")
    pdf.set_font("Helvetica", size=10)
    pdf.cell(0, 5, text, align="L", ln=1)
    return pdf


def WriteReport(
    workingdir: str,
    inpath: str,
    startpath: str,
    conf: configparser.ConfigParser,
    snakemake_resource_settings: ResourceSettings,
    snakemake_output_conf: OutputSettings,
    inputs_config: CLIparser,
    status: Literal["Failed", "Success"],
) -> None:
    if os.getcwd() != workingdir:
        os.chdir(workingdir)

    # sconfig.update(sparams)

    directories: dict[int, list[str]] = {
        0: [
            "Output directory:",
            "This is the directory where the output files were written as well as this summary.",
            "\t\t\t\t" + workingdir,
        ],
        1: [
            "Input directory:",
            "This is the directory that was given containing the input files.",
            "\t\t\t\t" + inpath,
        ],
        2: [
            "Starting point:",
            "This is the directory where the pipeline was started from.",
            "\t\t\t\t" + startpath,
        ],
    }

    pdf = PDF()
    pdf.set_author("SARS2seq")
    pdf.add_page()

    for i, k in enumerate(directories):
        pdf = directory_sections(pdf, i, directories)

    pdf.ln(10)

    if snakemake_output_conf.dryrun is True:
        pdf = analysis_details(pdf, "Workflow status:", "Dry-run")
    else:
        pdf = analysis_details(pdf, "Workflow status:", status)

    pdf.ln(5)

    # TODO: re-add the references, primers, and features sections to the report (if possible) when the rest of the pipeline compatible with the new methods

    # pdf = analysis_details(pdf, "Reference file:", sconfig["reference_file"])
    # pdf = analysis_details(pdf, "Primer file:", sconfig["primer_file"])
    # pdf = analysis_details(pdf, "GFF file:", sconfig["features_file"])
    pdf = analysis_details(pdf, "Sequencing platform:", inputs_config.flags.platform)
    pdf = analysis_details(
        pdf, "Selected amplicon type:", inputs_config.flags.amplicon_type
    )
    # if sconfig["primer_file"] != "None":
    #     pdf = analysis_details(
    #         pdf, "Primer mismatch rate:", str(sconfig["primer_mismatch_rate"])
    #     )

    pdf.ln(5)

    pdf = analysis_details(pdf, "Computing execution:", conf["COMPUTING"]["compmode"])
    if conf["COMPUTING"]["compmode"] == "grid":
        pdf = analysis_details(
            pdf, "Selected Grid Queue:", conf["COMPUTING"]["queuename"]
        )
    pdf = analysis_details(
        pdf, "Local available threads:", str(snakemake_resource_settings.cores)
    )

    pdf.ln(10)

    pdf = analysis_details(pdf, "Issued Command:", "")
    command = str(sys.argv[0]).split("/")[-1], *sys.argv[1:]
    pdf.multi_cell(0, 5, f'{" ".join(command)}')

    pdf.output(name="Runinfomation.pdf")

    log.debug("Python code :: Run report :: Create run report :: run information report created successfully.")
