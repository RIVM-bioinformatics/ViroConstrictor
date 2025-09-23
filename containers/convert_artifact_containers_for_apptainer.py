import json
import shutil
import subprocess
from typing import List

base_path_to_container_defs = "./containers"

if __name__ == "__main__":
    print("Renaming OCI containers and converting to Apptainer .sif format")

    builtcontainers: List = []
    with open(f"{base_path_to_container_defs}/builtcontainers.json", "r") as f:
        builtcontainers: List = json.load(f)

    builtcontainers_trimmed = [container.replace(":", "_") for container in builtcontainers]

    for original_name, trimmed_name in zip(builtcontainers, builtcontainers_trimmed):
        print(f"Renaming {original_name} to {trimmed_name}")
        shutil.move(
            f"{base_path_to_container_defs}/{original_name}.tar",
            f"{base_path_to_container_defs}/{trimmed_name}.tar",
        )
        trimmed_name_sif = trimmed_name.replace(":", "_")
        print(f"Converting {original_name} to Apptainer .sif format")
        subprocess.run(
            f"apptainer build {base_path_to_container_defs}/{trimmed_name_sif}.sif docker-archive://{base_path_to_container_defs}/{trimmed_name}.tar",
            shell=True,
        )
