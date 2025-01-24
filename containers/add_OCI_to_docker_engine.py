import json
import subprocess
from typing import List

base_path_to_container_defs = "./containers"

if __name__ == "__main__":
    print("Adding OCI containers to Docker Engine from Artifact")

    builtcontainers: List = []
    with open(f"{base_path_to_container_defs}/builtcontainers.json", "r") as f:
        builtcontainers: List = json.load(f)

    for container in builtcontainers:
        print(f"Adding {container} to local Docker Engine")
        subprocess.run(
            f"docker image load -i {base_path_to_container_defs}/{container}.tar",
            shell=True,
        )

    print("Done adding OCI containers to Docker Engine")
