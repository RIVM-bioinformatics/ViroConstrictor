import glob
import os
import pathlib
import shutil
import subprocess

from ViroConstrictor.userprofile import ReadConfig

userconf = ReadConfig(pathlib.Path("~/.ViroConstrictor_defaultprofile.ini").expanduser())

containerpath = None
try:
    containerpath = userconf["REPRODUCTION"]["container_cache_path"]
except KeyError as e:
    raise KeyError("The container_cache_path option is missing from the REPRODUCTION section of your configuration file.") from e

if not os.path.exists(containerpath):
    os.makedirs(containerpath)

subprocess.run(
    "python ./containers/build_containers.py",
    shell=True,
    check=True,
)

subprocess.run(
    "python containers/convert_artifact_containers_for_apptainer.py",
    shell=True,
    check=True,
)

for file in glob.glob("./containers/*.sif"):
    if os.path.exists(containerpath + "/" + file.split("/")[-1]):
        os.remove(containerpath + "/" + file.split("/")[-1])
    print(file)
    shutil.move(file, containerpath) if containerpath else None
for file in glob.glob("./containers/*.tar"):
    print(file)
    os.remove(file)

os.remove("./containers/builtcontainers.json")
