import json
import subprocess
from typing import List

# main_upstream_registry = "ghcr.io/rivm-bioinformatics"
main_upstream_registry = "ghcr.io/florianzwagemaker"

base_path_to_container_defs = "./containers"

if __name__ == "__main__":

    built_containers: List[str] = []
    with open(f"{base_path_to_container_defs}/builtcontainers.json", "r") as f:
        built_containers: List[str] = json.load(f)

    for container in built_containers:
        print(f"Tagging and pushing {container}")
        subprocess.run(
            f"docker tag {container} {main_upstream_registry}/{container}", shell=True
        )
        subprocess.run(f"docker push {main_upstream_registry}/{container}", shell=True)
