import os
import shutil
import sys

from ViroConstrictor.workflow.containers import download_containers

base_path_to_container_defs = "./containers"

if __name__ == "__main__":
    config = {"container_cache": f"{os.getcwd()}/test_containers", "dryrun": False}
    os.makedirs(config["container_cache"], exist_ok=True)

    # move .sif files from ./containers to ./test_containers
    for file in os.listdir(base_path_to_container_defs):
        if file.endswith(".sif"):
            shutil.move(
                f"{base_path_to_container_defs}/{file}",
                f"{config['container_cache']}/{file}",
            )

    download_status = download_containers(
        config["container_cache"], dryrun=config["dryrun"], verbose=True
    )
    if download_status == 1:
        print(
            "Failed to download all necessary containers. Please check the logs for more information."
        )
        sys.exit(1)
