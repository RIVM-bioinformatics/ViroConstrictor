import json
import os
import subprocess
import tempfile

import requests

from ViroConstrictor import __prog__
from ViroConstrictor.workflow.containers import fetch_hashes, upstream_registry

base_path_to_container_defs = "./containers"
upstream_api_endpoint = (
    "https://api.github.com/orgs/RIVM-bioinformatics/packages/container/"
)
upstream_api_authtoken = os.environ.get("TOKEN")
upstream_api_responsetype = "application/vnd.github+json"
upstream_api_version = "2022-11-28"
upstream_api_headers = {
    "Accept": f"{upstream_api_responsetype}",
    "X-GitHub-Api-Version": f"{upstream_api_version}",
    "Authorization": f"Bearer {upstream_api_authtoken}",
}


#TODO: break up this script into smaller functions
if __name__ == "__main__":
    print("Start of container building process for ViroConstrictor")
    recipe_hashes = fetch_hashes()

    builtcontainers = []
    for recipe, VersionHash in recipe_hashes.items():
        # strip the name of the recipe to only get the name of the environment
        recipe_basename = os.path.basename(recipe).replace(".yaml", "")
        container_basename = f"{__prog__}_{recipe_basename}".lower()

        associated_container_dock_file = os.path.join(
            base_path_to_container_defs, f"{recipe_basename}.dockerfile"
        )
        upstream_registry_url = f"{upstream_registry}/{recipe_basename}:{VersionHash}"
        upstream_existing_containers = (
            f"{upstream_api_endpoint}{__prog__}_{recipe_basename}/versions"
        )
        print(
            f"Checking if container '{container_basename}' with hash '{VersionHash}' exists in the upstream registry"
        )
        json_response = requests.get(
            upstream_existing_containers, headers=upstream_api_headers
        ).json()

        tags = []

        # if the container exists at all in the upstream registry, the json response will be a list.
        # If the container does not exist, the json response will be a dict with a message that the container does not exist.
        # You can therefore check if the json response is a list or a dict to see if the container exists or not.
        if isinstance(json_response, list):
            tags = [
                version["metadata"]["container"]["tags"] for version in json_response
            ]
            # flatten the list of tags
            tags = [tag for sublist in tags for tag in sublist]
            print(tags)

        if VersionHash in tags:
            print(
                f"Container '{container_basename}' with hash '{VersionHash}' already exists in the upstream registry"
            )
            continue

        print(
            f"Container '{container_basename}' with hash '{VersionHash}' does not exist in the upstream registry"
        )
        print(
            f"Starting Docker build process for container '{container_basename}:{VersionHash}'"
        )

        # create a temporary file to write the container definition to, copy the contents of {recipe_basename}.dockerfile to it and then append the labels section to it including the version hash
        # then use the temporary file as the container definition file for the docker build process
        # the docker build process will build the container file also in a temporary directory
        # after the container is built, the built container will be saved in the docker artifact database (local).
        # This is necessary to transform the container from docker format to apptainer format in a separate script.
        # the container file will not be pushed to the upstream registry yet, this will be done in a separate script after all containers have been built and tested.
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False
        ) as tmp, tempfile.TemporaryDirectory() as tmpdir:
            with open(associated_container_dock_file, "r") as f:
                tmp.write(f.read())
                tmp.write(
                    f"""

LABEL Author="RIVM-bioinformatics team"
LABEL Maintainer="RIVM-bioinformatics team"
LABEL Associated_pipeline="{__prog__}"
LABEL version="{VersionHash}"
LABEL org.opencontainers.image.authors="ids-bioinformatics@rivm.nl"
LABEL org.opencontainers.image.source=https://github.com/RIVM-bioinformatics/{__prog__}

    """
                )
            tmp.flush()  # flush the temporary file to make sure the contents are written to disk
            subprocess.run(
                [
                    "docker",
                    "build",
                    "-t",
                    f"{container_basename}:{VersionHash}",
                    "-f",
                    f"{tmp.name}",
                    ".",
                    "--network",
                    "host",
                    "--no-cache",
                ],
                check=True,
            )
            # move the container file to the current working directory
            subprocess.run(
                [
                    "docker",
                    "save",
                    "-o",
                    f"{base_path_to_container_defs}/{container_basename}:{VersionHash}.tar",
                    f"{container_basename}:{VersionHash}",
                ]
            )

        builtcontainers.append(f"{container_basename}:{VersionHash}")

    with open(f"{base_path_to_container_defs}/builtcontainers.json", "w") as f:
        json.dump(builtcontainers, f, indent=4)
