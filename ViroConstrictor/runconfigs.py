# pylint: disable=C0103

"""
Construct and write configuration files for SnakeMake
"""

import multiprocessing
import os

import yaml


def WriteYaml(data, filepath):
    '''WriteYaml takes a dictionary and a filepath, and writes the given dictionary to the filepath as a yaml
    file
    
    Parameters
    ----------
    data
        The data to be written to the file.
    filepath
        The path to the file you want to write to.
    
    Returns
    -------
        The filepath
    
    '''
    if not os.path.exists(os.path.dirname(filepath)):
        os.makedirs(os.path.dirname(filepath))
    with open(filepath, "w") as file:
        yaml.dump(data, file, default_flow_style=False)
    return filepath


def set_cores(cores):
    '''Balance the requested number of cores to the available number of cores on the system to ensure analysis will not overload the system in local-analysis mode.
    If the number of cores requested is equal to the number of cores available, return the number of cores requested minus 2.
    If the number of cores requested is greater than the number of cores available, then return the
    number of cores available minus 2. 
    Otherwise, return the number of cores requested
    
    Parameters
    ----------
    cores
        The number of cores to use for parallel processing given by the user.
    
    Returns
    -------
        The number of cores that can be used by ViroConstrictor during analysis.
    
    '''
    available = multiprocessing.cpu_count()
    if cores == available:
        return cores - 2
    if cores > available:
        return available - 2
    return cores


def get_max_local_mem():
    '''It returns the maximum amount of local memory that can be allocated to a single process on the
    current machine
    
    Returns
    -------
        The amount of memory available in MB.
    
    '''
    avl_mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    return int(round(avl_mem_bytes / (1024.0 ** 2) - 2000, -3))


def SnakemakeConfig(conf, cores, dryrun):
    '''Function takes the configuration file and the number of cores to use, and returns a dictionary that can be
    used to configure Snakemake
    
    Parameters
    ----------
    conf
        the configuration file
    cores
        The number of cores to use for each job.
    dryrun
        If True, the pipeline will only execute in dry-run mode.
    
    Returns
    -------
        A dictionary with the configuration for the snakemake workflow.
    
    '''
    cores = set_cores(cores)
    compmode = conf["COMPUTING"]["compmode"]

    if compmode == "local":
        config = {
            "cores": cores,
            "latency-wait": 60,
            "use-conda": True,
            "dryrun": dryrun,
            "jobname": "ViroConstrictor_{name}.jobid{jobid}",
        }

    if compmode == "grid":
        queuename = conf["COMPUTING"]["queuename"]
        threads = "{threads}"
        mem = "{resources.mem_mb}"
        config = {
            "cores": 300,
            "latency-wait": 60,
            "use-conda": True,
            "dryrun": dryrun,
            "jobname": "ViroConstrictor_{name}.jobid{jobid}",
            "drmaa": f' -q {queuename} -n {threads} -R "span[hosts=1]" -M {mem}',
            "drmaa-log-dir": "logs/drmaa",
        }

    return config

# Todo: possibly refactor this function to reduce parameters
def SnakemakeParams(conf, cores, sampleinfo, platform, samplesheet, amplicon_type):
    '''Function takes the configuration file, the number of cores, the sample information, the platform, the
    sample sheet, and the amplicon type, and returns a dictionary with the sample sheet, the computing
    execution mode, the maximum local memory, the platform, the amplicon type, the number of threads for
    each step, and the run parameters
    
    Parameters
    ----------
    conf
        the configuration file
    cores
        the number of cores to use for the pipeline
    sampleinfo
        a dictionary of sample information
    platform
        The sequencing platform used to generate the data.
    samplesheet
        the path to the samplesheet
    amplicon_type
        This is the type of amplicon that you are using. This is used to determine the primer sequences to
    remove.
    
    Returns
    -------
        A dictionary with the following keys:
        sample_sheet
        computing_execution
        max_local_mem
        platform
        amplicon_type
        threads
        runparams
    
    '''
    if conf["COMPUTING"]["compmode"] == "local":
        threads_highcpu = int(cores - 2)
        threads_midcpu = int(cores / 2)
        threads_lowcpu = 1
    if conf["COMPUTING"]["compmode"] == "grid":
        threads_highcpu = 12
        threads_midcpu = 6
        threads_lowcpu = 2

    return {
        "sample_sheet": samplesheet,
        "computing_execution": conf["COMPUTING"]["compmode"],
        "max_local_mem": get_max_local_mem(),
        "platform": platform,
        "amplicon_type": amplicon_type,
        "threads": {
            "Alignments": threads_highcpu,
            "QC": threads_midcpu,
            "AdapterRemoval": threads_lowcpu,
            "PrimerRemoval": threads_highcpu,
            "Consensus": threads_midcpu,
            "Index": threads_lowcpu,
            "Typing": threads_lowcpu,
        },
        "runparams": {
            "alignmentfilters": "-F 256 -F 512 -F 4 -F 2048",
            "qc_filter_illumina": 20,
            "qc_filter_nanopore": 7,
            "qc_filter_iontorrent": 20,
            "qc_window_illumina": 5,
            "qc_window_nanopore": 20,
            "qc_window_iontorrent": 15,
            "qc_min_readlength": 100,
        },
    }


def LoadConf(configfile):
    '''Function opens a yaml file, reads it, and returns the content as a dictionary
    
    Parameters
    ----------
    configfile
        The path to the yaml file.
    
    Returns
    -------
        A dictionary of the yaml file.
    
    '''
    with open(configfile, "r") as ConfIn:
        conf = yaml.load(ConfIn, Loader=yaml.FullLoader)
    return conf
