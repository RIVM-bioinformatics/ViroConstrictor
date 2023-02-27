# pylint: disable=C0103
"""
Write the samplesheets
"""

import os
import pathlib
import re


def illumina_sheet(inputdir: pathlib.Path) -> dict[str, dict[str, str]]:
    """Takes a directory, walks through it, and returns a dictionary of dictionaries

    Parameters
    ----------
    inputdir : pathlib.Path
        pathlib.Path = pathlib.Path("/path/to/input/directory")

    Returns
    -------
        A dictionary of dictionaries. The outer dictionary has the sample name as the key and the inner
    dictionary has the read number as the key and the file path as the value.

    """
    illuminapattern: re.Pattern = re.compile(
        r"(.*)(_|\.)R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?"
    )
    samples: dict[str, dict[str, str]] = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath: str = os.path.abspath(os.path.join(dirname, files))
            if match := illuminapattern.fullmatch(files):
                sample = samples.setdefault(match[1], {})
                sample[f"R{match[3]}"] = fullpath
    return samples


def nanopore_sheet(inputdir: pathlib.Path) -> dict[str, str]:
    """Takes a directory as input, and returns a dictionary of sample names and their corresponding file
    paths

    Parameters
    ----------
    inputdir : pathlib.Path
        pathlib.Path = pathlib.Path("/path/to/input/directory")

    Returns
    -------
        A dictionary of samples and their full path.

    """
    nanoporepattern: re.Pattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples: dict[str, str] = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath: str = os.path.abspath(os.path.join(dirname, files))
            if match := nanoporepattern.fullmatch(files):
                samples.setdefault(match[1], fullpath)
    return samples


def iontorrent_sheet(inputdir: pathlib.Path) -> dict[str, str]:
    """Takes a directory as input, and returns a dictionary of sample names and their corresponding file
    paths

    Parameters
    ----------
    inputdir : pathlib.Path
        pathlib.Path = pathlib.Path("/path/to/input/directory")

    Returns
    -------
        A dictionary with the sample name as the key and the full path to the file as the value.

    """
    iontorrentpattern: re.Pattern = re.compile(r"(.*)\.f(ast)?q(\.gz)?")
    samples: dict[str, str] = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath: str = os.path.abspath(os.path.join(dirname, files))
            if match := iontorrentpattern.fullmatch(files):
                samples.setdefault(match[1], fullpath)
    return samples


def GetSamples(
    inputdir: pathlib.Path, platform: str
) -> dict[str, str] | dict[str, dict[str, str]]:
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
        illumina_samples: dict[str, dict[str, str]] = illumina_sheet(inputdir)
        return illumina_samples
    elif platform == "iontorrent":
        iontorrent_samples: dict[str, str] = iontorrent_sheet(inputdir)
        return iontorrent_samples
    elif platform == "nanopore":
        nanopore_samples: dict[str, str] = nanopore_sheet(inputdir)
        return nanopore_samples
    return {}
