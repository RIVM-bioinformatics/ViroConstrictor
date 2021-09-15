# pylint: disable=C0103
"""
Write the samplesheets
"""

import os
import re

import yaml


def illumina_sheet(inputdir, sheet):
    illuminapattern = re.compile(r"(.*)(_|\.)R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = illuminapattern.fullmatch(files)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["R{}".format(match.group(3))] = str(fullpath)
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()


def nanopore_sheet(inputdir, sheet):
    nanoporepattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = nanoporepattern.fullmatch(files)
            if match:
                samples.setdefault(match.group(1), fullpath)
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()


def iontorrent_sheet(inputdir, sheet):
    nanoporepattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = nanoporepattern.fullmatch(files)
            if match:
                samples.setdefault(match.group(1), fullpath)
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()


def WriteSampleSheet(inputdir, platform):
    if platform == "illumina":
        illumina_sheet(inputdir, "samplesheet.yaml")
    if platform == "nanopore":
        nanopore_sheet(inputdir, "samplesheet.yaml")
    if platform == "iontorrent":
        iontorrent_sheet(inputdir, "samplesheet.yaml")

    samplesheet = os.getcwd() + "/samplesheet.yaml"
    return samplesheet
