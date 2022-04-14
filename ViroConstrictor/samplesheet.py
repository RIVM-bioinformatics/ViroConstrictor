# pylint: disable=C0103
"""
Write the samplesheets
"""

import os
import re

import yaml


def illumina_sheet(inputdir):
    """Function takes a directory as input, and returns a dictionary of dictionaries, where the keys of the outer
    dictionary are the sample names, and the keys of the inner dictionaries are the read numbers (R1 and
    R2)
    
    Parameters
    ----------
    inputdir
        The directory where the fastq files are located.
    
    Returns
    -------
        A dictionary of dictionaries.
    
    """
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
    """Function takes a directory as input, and returns a dictionary of sample names and their corresponding file
    paths
    
    Parameters
    ----------
    inputdir
        the directory where the fastq files are located
    
    Returns
    -------
        A dictionary with the sample name as the key and the full path to the file as the value.
    
    """
    nanoporepattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.abspath(os.path.join(dirname, files))
            if match := nanoporepattern.fullmatch(files):
                samples.setdefault(match[1], fullpath)
    return samples


def iontorrent_sheet(inputdir):
    """Function takes a directory as input, and returns a dictionary of sample names and their corresponding file
    paths
    
    Parameters
    ----------
    inputdir
        the directory where the fastq files are located
    
    Returns
    -------
        A dictionary with the sample name as the key and the full path to the file as the value.
    
    """
    iontorrentpattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.abspath(os.path.join(dirname, files))
            if match := iontorrentpattern.fullmatch(files):
                samples.setdefault(match[1], fullpath)
    return samples


def GetSamples(inputdir, platform):
    """Wrapping function taking in a directory and sequencing platform, triggers appropriate sub-function and returns a dictionary of samples
    
    Parameters
    ----------
    inputdir
        the directory where the sample sheets are located
    platform
        the sequencing platform used to generate the data.
    
    Returns
    -------
        A dict of samples
    
    """
    if platform == "illumina":
        samples = illumina_sheet(inputdir)
    if platform == "nanopore":
        samples = nanopore_sheet(inputdir)
    if platform == "iontorrent":
        samples = iontorrent_sheet(inputdir)
    return samples
