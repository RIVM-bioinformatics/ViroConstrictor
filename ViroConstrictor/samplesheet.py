# pylint: disable=C0103
"""
Write the samplesheets
"""

import os
import re

import yaml


def illumina_sheet(inputdir):
    illuminapattern = re.compile(r"(.*)(_|\.)R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = illuminapattern.fullmatch(files)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["R{}".format(match.group(3))] = str(fullpath)
    return samples
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()


def nanopore_sheet(inputdir):
    nanoporepattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = nanoporepattern.fullmatch(files)
            if match:
                samples.setdefault(match.group(1), fullpath)
    return samples
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()


def iontorrent_sheet(inputdir):
    iontorrentpattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = iontorrentpattern.fullmatch(files)
            if match:
                samples.setdefault(match.group(1), fullpath)
    return samples
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()

def GetSamples(inputdir, platform):
    if platform == "illumina":
        samples = illumina_sheet(inputdir)
    if platform == "nanopore":
        samples = nanopore_sheet(inputdir)
    if platform == "iontorrent":
        samples = iontorrent_sheet(inputdir)
    return samples

def WriteSampleSheet(inputdir, platform):
    samples = GetSamples(inputdir, platform)
    with open("samplesheet.yaml", "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)

    samplesheet = f'{os.getcwd()}/samplesheet.yaml'
    return samplesheet
