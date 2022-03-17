import os
import sys
from datetime import datetime

from fpdf import FPDF

from ViroConstrictor import __version__


class PDF(FPDF):
    def timestamp(self):
        return datetime.now().strftime("%d-%m-%Y %H:%M")

    def header(self):
        self.set_font("Helvetica", size=23, style="B")
        self.cell(0, 20, "ViroConstrictor", align="C", ln=1)
        self.set_font("Arial", size=12)
        self.cell(0, 5, "Run information summary", align="C", ln=1)
        self.set_font("Arial", size=8)
        self.cell(0, 5, self.timestamp(), align="C", ln=1)
        self.set_font("Arial", size=10, style="BU")
        self.cell(0, 15, "Version: " + __version__, align="C", ln=1)


def directory_sections(pdf, iteration, contents):
    pdf.set_font("Arial", size=12, style="B")
    pdf.cell(40, 12, contents.get(iteration)[0], align="L")
    pdf.set_font("Arial", size=10)
    pdf.cell(0, 12, contents.get(iteration)[1], align="L", ln=1)
    pdf.set_font("Arial", size=12, style="I")
    pdf.cell(0, 5, contents.get(iteration)[2], align="L", ln=1)
    return pdf


def analysis_details(pdf, header, text):
    pdf.set_font("Arial", size=12, style="B")
    pdf.cell(55, 5, header, align="L")
    pdf.set_font("Arial", size=10)
    pdf.cell(0, 5, text, align="L", ln=1)
    return pdf


def WriteReport(workingdir, inpath, startpath, conf, sparams, sconfig, status):
    if os.getcwd() != workingdir:
        os.chdir(workingdir)

    sconfig.update(sparams)

    directories = {
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

    pdf = analysis_details(pdf, "Reference file:", sconfig["reference_file"])
    pdf = analysis_details(pdf, "Primer file:", sconfig["primer_file"])
    pdf = analysis_details(pdf, "GFF file:", sconfig["features_file"])
    pdf = analysis_details(pdf, "Sequencing platform:", sconfig["platform"])
    pdf = analysis_details(pdf, "Selected amplicon type:", sconfig["amplicon_type"])
    if sconfig["primer_file"] != "None":
        pdf = analysis_details(
            pdf, "Primer mismatch rate:", str(sconfig["primer_mismatch_rate"])
        )

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

    pdf.output("Runinfomation.pdf", "F")
