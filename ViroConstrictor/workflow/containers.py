import hashlib
import os
import subprocess
from typing import Any, Dict, List, Tuple

from rich import print

from ViroConstrictor import __prog__

upstream_registry = "ghcr.io/rivm-bioinformatics"


def fetch_recipes(recipe_folder: str) -> List[str]:
    return [
        os.path.abspath(os.path.join(recipe_folder, file))
        for file in os.listdir(recipe_folder)
        if file.endswith(".yaml")
    ]


def fetch_scripts(script_folder: str) -> List[str]:
    script_files = []
    for root, dirs, files in os.walk(script_folder):
        script_files.extend(
            os.path.abspath(os.path.join(root, file))
            for file in files
            if file.endswith(".py")
        )
    return script_files


def fetch_files(file_folder: str) -> List[str]:
    return [
        os.path.abspath(os.path.join(file_folder, file))
        for file in os.listdir(file_folder)
    ]


def fetch_hashes() -> dict[str, str]:
    recipe_files = sorted(
        fetch_recipes(f"{os.path.dirname(os.path.realpath(__file__))}/envs/")
    )
    script_files = sorted(
        fetch_scripts(f"{os.path.dirname(os.path.realpath(__file__))}/scripts/")
    )
    config_files = sorted(
        fetch_files(f"{os.path.dirname(os.path.realpath(__file__))}/files/")
    )

    script_hashes = {}
    for script_file in script_files:
        with open(script_file, "rb") as f:
            file_hash = hashlib.sha256(f.read()).hexdigest()[:6]
            script_hashes[script_file] = file_hash

    config_hashes = {}
    for config_file in config_files:
        with open(config_file, "rb") as f:
            file_hash = hashlib.sha256(f.read()).hexdigest()[:6]
            config_hashes[config_file] = file_hash

    # sort the hashes of the scripts and the configs
    script_hashes = dict(sorted(script_hashes.items()))
    config_hashes = dict(sorted(config_hashes.items()))

    # join the hashes of the scripts and the configs (the values of the dictionaries), make a new hash of the joined hashes
    merged_hashes = hashlib.sha256(
        "".join(list(script_hashes.values()) + list(config_hashes.values())).encode()
    ).hexdigest()[:6]

    hashes = {}
    for recipe_file in recipe_files:
        with open(recipe_file, "rb") as f:
            recipe_hash = hashlib.sha256(f.read()).hexdigest()[:6]
            if os.path.basename(recipe_file).split(".")[0] != "Scripts":
                hashes[recipe_file] = recipe_hash
                continue
            # add the merged hash to the recipe hash and make a new hash of the joined hashes
            file_hash = hashlib.sha256(
                (recipe_hash + merged_hashes).encode()
            ).hexdigest()[:6]
            hashes[recipe_file] = file_hash
    return hashes


def get_hash(target_container: str) -> str | None:
    hashes = fetch_hashes()
    return next(
        (hash for recipe, hash in hashes.items() if target_container in recipe),
        None,
    )


def containerization_installed() -> bool:
    if os.system("which apptainer") == 0:
        return True
    return os.system("which singularity") == 0


def containers_to_download(config: Dict[str, Any]) -> List[str]:
    recipes = fetch_hashes()

    # check the recipes dict and create a list of required containers for the workflow. the list must contain strings that look like "viroconstrictor_alignment_{hash}" "viroconstrictor_clean_{hash}" etc.
    required_containers = []
    for key, val in recipes.items():
        recipe_basename = os.path.basename(key).replace(".yaml", "")
        required_containers.append(f"{__prog__}_{recipe_basename}_{val}".lower())

    # check if the required containers are already present in the directory listed in config["container_cache"]
    # if they are, remove them from the list
    containers_present = os.listdir(config["container_cache"])

    # remove the .sif extension from the container names
    containers_present = [item.split(".")[0] for item in containers_present]

    containers_to_download = []
    for container in required_containers:
        if container in containers_present:
            continue
        containers_to_download.append(container)
    return containers_to_download


def containerization_executable() -> str:
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


def download_containers(config: Dict[str, Any]) -> int:
    to_download = containers_to_download(config)
    to_download = [x.rsplit("_", 1)[0] + ":" + x.rsplit("_", 1)[1] for x in to_download]

    # I thought this would be a great place to use concurrent.futures to download the containers in parallel
    # however this has unintended consequences for the various OCI layers resulting in at least one container not being downloaded correctly.
    # Thus it's better for now to download the containers sequentially.
    for container in to_download:
        print(f"Downloading {container}")
        executable = containerization_executable()
        status = subprocess.call(
            f"{executable} pull --dir {config['container_cache']} docker://{upstream_registry}/{container}",
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        if status != 0:
            print(f"Failed to download {container}")
            return 1
        print(f"Successfully downloaded {container}")

    return 0


def construct_container_bind_args(samples_dict: Dict) -> str:
    paths = []
    for nested_dict in samples_dict.items():
        paths.extend(
            f"{os.path.dirname(value)}"
            for value in nested_dict.values()
            if isinstance(value, str) and os.path.exists(value)
        )
    # remove all duplicates from the paths list by making it a set
    # for every item in the set, add '--bind '
    # join all the items together to make a long string
    return " ".join([f"--bind {path}" for path in set(paths)])
