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
            fullpath = os.path.abspath(os.path.join(dirname, files))
            if match := illuminapattern.fullmatch(files):
                sample = samples.setdefault(match[1], {})
                sample[f"R{match[3]}"] = str(fullpath)
    return samples


def nanopore_sheet(inputdir):
    nanoporepattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.abspath(os.path.join(dirname, files))
            if match := nanoporepattern.fullmatch(files):
                samples.setdefault(match[1], fullpath)
    return samples


def iontorrent_sheet(inputdir):
    iontorrentpattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.abspath(os.path.join(dirname, files))
            if match := iontorrentpattern.fullmatch(files):
                samples.setdefault(match[1], fullpath)
    return samples


def GetSamples(inputdir, platform):
    if platform == "illumina":
        samples = illumina_sheet(inputdir)
    if platform == "nanopore":
        samples = nanopore_sheet(inputdir)
    if platform == "iontorrent":
        samples = iontorrent_sheet(inputdir)
    return samples
