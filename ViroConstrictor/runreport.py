import configparser
import os
import sys
from datetime import datetime
from typing import Any, Literal

from fpdf import FPDF

from ViroConstrictor import __version__


class PDF(FPDF):
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
    sparams: dict[str, Any],
    sconfig: dict[str, Any],
    status: Literal["Failed", "Success"],
) -> None:
    if os.getcwd() != workingdir:
        os.chdir(workingdir)

    sconfig.update(sparams)

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

    if sconfig["dryrun"] is True:
        pdf = analysis_details(pdf, "Workflow status:", "Dry-run")
    else:
        pdf = analysis_details(pdf, "Workflow status:", status)

    pdf.ln(5)

    # TODO: re-add the references, primers, and features sections to the report (if possible) when the rest of the pipeline compatible with the new methods

    # pdf = analysis_details(pdf, "Reference file:", sconfig["reference_file"])
    # pdf = analysis_details(pdf, "Primer file:", sconfig["primer_file"])
    # pdf = analysis_details(pdf, "GFF file:", sconfig["features_file"])
    pdf = analysis_details(pdf, "Sequencing platform:", sconfig["platform"])
    pdf = analysis_details(pdf, "Selected amplicon type:", sconfig["amplicon_type"])
    # if sconfig["primer_file"] != "None":
    #     pdf = analysis_details(
    #         pdf, "Primer mismatch rate:", str(sconfig["primer_mismatch_rate"])
    #     )

    pdf.ln(5)

    pdf = analysis_details(pdf, "Computing execution:", sconfig["computing_execution"])
    if sconfig["computing_execution"] == "grid":
        pdf = analysis_details(
            pdf, "Selected Grid Queue:", conf["COMPUTING"]["queuename"]
        )
    pdf = analysis_details(pdf, "Local available threads:", str(sconfig["cores"]))

    pdf.ln(10)

    pdf = analysis_details(pdf, "Issued Command:", "")
    command = str(sys.argv[0]).split("/")[-1], *sys.argv[1:]
    pdf.multi_cell(0, 5, f'{" ".join(command)}')

    pdf.output(name="Runinfomation.pdf")
