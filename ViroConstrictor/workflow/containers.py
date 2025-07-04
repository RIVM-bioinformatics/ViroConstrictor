import hashlib
import os
import subprocess
from pathlib import Path
from typing import Dict, List

from ViroConstrictor import __prog__
from ViroConstrictor.logging import log

upstream_registry = "ghcr.io/rivm-bioinformatics"


def fetch_recipes(recipe_folder: str) -> List[str]:
    """
    Fetches a list of absolute file paths for all YAML files in the given recipe folder.

    Parameters
    ----------
    recipe_folder : str
        The path to the folder containing the recipe files.

    Returns
    -------
    List[str]
        A list of absolute file paths for all YAML files in the recipe folder.
    """
    return [
        os.path.abspath(os.path.join(recipe_folder, file))
        for file in os.listdir(recipe_folder)
        if file.endswith(".yaml")
    ]


def fetch_scripts(script_folder: str) -> List[str]:
    """
    Fetches all Python script files with the extension '.py' from the specified folder.

    Parameters
    ----------
    script_folder : str
        The path to the folder containing the script files.

    Returns
    -------
    List[str]
        A list of absolute paths to the fetched script files.
    """
    script_files = []
    for root, dirs, files in os.walk(script_folder):
        script_files.extend(
            os.path.abspath(os.path.join(root, file))
            for file in files
            if file.endswith(".py")
        )
    return script_files


def fetch_files(file_folder: str) -> List[str]:
    """
    Fetches a list of absolute file paths from the given file folder.

    Parameters
    ----------
    file_folder : str
        The path to the folder containing the files.

    Returns
    -------
    List[str]
        A list of absolute file paths.

    """
    return [
        os.path.abspath(os.path.join(file_folder, file))
        for file in os.listdir(file_folder)
    ]


def calculate_hashes(file_list: List[str]) -> Dict[str, str]:
    """
    Parameters
    ----------
    file_list : List[str]
        A list of file paths for which to calculate the hashes.

    Returns
    -------
    Dict[str, str]
        A dictionary where the keys are file paths and the values are the first 6 characters of the SHA-256 hash of the file contents.
    """
    hashdict = {}
    for file in file_list:
        with open(file, "rb") as f:
            file_hash = hashlib.sha256(f.read()).hexdigest()[:6]
            hashdict[file] = file_hash
    return hashdict


def fetch_hashes() -> Dict[str, str]:
    """
    Fetches and returns the hashes of recipe files, script files, and config files.

    Returns
    -------
    Dict[str, str]
        A dictionary containing the file paths as keys and their corresponding hashes as values.
    """
    # Fetch the recipe files, script files, and config files
    recipe_files = sorted(
        fetch_recipes(f"{os.path.dirname(os.path.realpath(__file__))}/envs/")
    )
    script_files = sorted(
        fetch_scripts(f"{os.path.dirname(os.path.realpath(__file__))}/scripts/")
    )
    config_files = sorted(
        fetch_files(f"{os.path.dirname(os.path.realpath(__file__))}/files/")
    )

    # Calculate hashes for script files
    script_hashes = calculate_hashes(script_files)

    # Calculate hashes for config files
    config_hashes = calculate_hashes(config_files)

    # Sort the hashes of the scripts and the configs
    script_hashes = dict(sorted(script_hashes.items()))
    config_hashes = dict(sorted(config_hashes.items()))

    # Join the hashes of the scripts and the configs, and create a new hash of the joined hashes
    merged_hashes = hashlib.sha256(
        "".join(list(script_hashes.values()) + list(config_hashes.values())).encode()
    ).hexdigest()[:6]

    # Calculate hashes for recipe files
    hashes = {}
    for recipe_file in recipe_files:
        with open(recipe_file, "rb") as f:
            recipe_hash = hashlib.sha256(f.read()).hexdigest()[:6]
            # if the recipe file is not named 'Scripts', then add the hash to the dictionary and continue
            if os.path.basename(recipe_file).split(".")[0] != "Scripts":
                hashes[recipe_file] = recipe_hash
                continue
            # if the recipe file *is* named 'Scripts', then combine the recipe hash with the merged hash of the scripts and configs
            # Subsequently, create a new hash of the joined hashes to avoid conflicts and end up with a singular tracker.
            file_hash = hashlib.sha256(
                (recipe_hash + merged_hashes).encode()
            ).hexdigest()[:6]
            hashes[recipe_file] = file_hash

    return hashes


def get_hash(target_container: str) -> str | None:
    """
    Get the hash for the target container.

    Parameters
    ----------
    target_container : str
        The name of the target container.

    Returns
    -------
    str or None
        The hash of the target container if found, None otherwise.
    """
    hashes = fetch_hashes()
    return next(
        (hash for recipe, hash in hashes.items() if target_container in recipe),
        None,
    )


def containerization_installed() -> bool:
    """
    Check if containerization tools are installed.

    Returns
    -------
    bool
        True if either 'apptainer' or 'singularity' is installed, False otherwise.
    """
    if os.system("which apptainer") == 0:
        return True
    return os.system("which singularity") == 0


def containers_to_download(apptainer_path: Path) -> List[str]:
    """
    Returns a list of containers that need to be downloaded for the workflow.

    Parameters
    ----------
    apptainer_path : Path
        The path to the Apptainer container directory.

    Returns
    -------
    List[str]
        A list of container names that need to be downloaded.

    """
    recipes = fetch_hashes()

    # check the recipes dict and create a list of required containers for the workflow. the list must contain strings that look like "viroconstrictor_alignment_{hash}" "viroconstrictor_clean_{hash}" etc.
    required_containers = []
    for key, val in recipes.items():
        recipe_basename = os.path.basename(key).replace(".yaml", "")
        required_containers.append(f"{__prog__}_{recipe_basename}_{val}".lower())

    # check if the folder exists, if not create it
    if apptainer_path is not None and not os.path.exists(apptainer_path):
        os.makedirs(apptainer_path, exist_ok=True)
    containers_present = os.listdir(apptainer_path)

    # remove the .sif extension from the container names
    containers_present = [item.split(".")[0] for item in containers_present]

    # loop through the required_containers list and check if they are present in the containers_present list
    # if they are not present, add them to the containers_to_download list
    return [
        container
        for container in required_containers
        if container not in containers_present
    ]


def containerization_executable() -> str:
    """
    Determines the containerization executable to be used.

    Returns
    -------
    str
        The name of the containerization executable ('apptainer' or 'singularity').
    """
    return (
        "apptainer"
        if subprocess.call(
            "which apptainer",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        == 0
        else "singularity"
    )


def download_containers(
    apptainer_path: Path, dryrun: bool = False, verbose=False
) -> int:
    """
    Download containers specified in the configuration.

    Parameters
    ----------
    apptainer_path : Path
        The path to the Apptainer container directory.
    dryrun : bool, optional
        If True, only simulate the download without actually performing it. Defaults to False.
    verbose : bool, optional
        Whether to display verbose output. Defaults to False.

    Returns
    -------
    int
        0 if all containers were downloaded successfully, 1 otherwise.
    """
    to_download = containers_to_download(apptainer_path)
    to_download = [x.rsplit("_", 1)[0] + ":" + x.rsplit("_", 1)[1] for x in to_download]

    if dryrun:
        log.info(
            f"Container(s) [magenta]{', '.join(to_download)}[/magenta] will be downloaded"
        )
        return 0
    # I thought this would be a great place to use concurrent.futures to download the containers in parallel
    # however this has unintended consequences for the various OCI layers resulting in at least one container not being downloaded correctly.
    # Thus it's better for now to download the containers sequentially.
    for container in to_download:
        log.info(
            f"Downloading container: [magenta]'{container}'[/magenta] to local cache"
        )
        executable = containerization_executable()
        status = subprocess.call(
            f"{executable} pull --dir {apptainer_path} docker://{upstream_registry}/{container}",
            shell=True,
            stderr=subprocess.PIPE if verbose is False else None,
            stdout=subprocess.PIPE if verbose is False else None,
        )
        if status != 0:
            log.error(f"Failed to download container: [magenta]'{container}'[/magenta]")
            return 1
        log.info(
            f"Successfully downloaded container: [magenta]'{container}'[/magenta] to local cache"
        )

    return 0


def construct_container_bind_args(samples_dict: Dict) -> str:
    """
    Constructs the bind arguments for container execution based on the given samples dictionary.

    Parameters
    ----------
    samples_dict : Dict
        A dictionary containing sample information.

    Returns
    -------
    str
        A string representing the bind arguments for container execution.

    Notes
    -----
    This function iterates over the given samples dictionary and extracts the paths of all files that exist.
    It then removes any duplicate paths and constructs a string of bind arguments for container execution.

    Examples
    --------
    >>> samples = {
    ...     'sample1': {
    ...         'file1': '/path/to/file1',
    ...         'file2': '/path/to/file2'
    ...     },
    ...     'sample2': {
    ...         'file3': '/path/to/file3',
    ...         'file4': '/path/to/file4'
    ...     }
    ... }
    >>> construct_container_bind_args(samples)
    '--bind /path/to'
    """
    paths = []
    for keys, nested_dict in samples_dict.items():
        paths.extend(
            f"{os.path.dirname(value)}"
            for value in nested_dict.values()
            if isinstance(value, str) and os.path.exists(value)
        )
    # remove all duplicates from the paths list by making it a set
    # for every item in the set, add '--bind '
    # join all the items together to make a long string
    return " ".join([f"--bind {path}" for path in set(paths)])
