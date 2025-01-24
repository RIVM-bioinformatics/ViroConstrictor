# Installation

ViroConstrictor is only available on Linux (or Linux-based) operating systems. MacOS may also work but is not tested.
ViroConstrictor will *not work* on Windows.  


## Prerequisites

ViroConstrictor requires the following software to be installed on your system:

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://mamba.readthedocs.io/en/latest/installation.html)

### Recommended but optional

- [Apptainer](https://apptainer.org/)


## Installation with conda

ViroConstrictor releases are published on Bioconda to allow for easy installation.

If you already have [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) installed then please use the command below:
```bash
mamba install -c conda-forge -c bioconda viroconstrictor
```

??? info "Click here if you don't have mamba installed in your system already"
    ```bash
    conda install -c conda-forge -c bioconda viroconstrictor
    ```

Installation through mamba/conda is the recommended method for installing ViroConstrictor.

!!! tip "Use a dedicated conda environment for ViroConstrictor"
    It is recommended to install and use ViroConstrictor in a dedicated conda-environment. This is to ensure there will be no conflicting dependencies.

    You can create a dedicated environment and install viroconstrictor in it through a single command with both mamba and conda.

    The following command will use **mamba** to create a new environment named 'pipeline_env' with viroconstrictor installed in it.
    ```bash
    mamba create --name pipeline_env -c conda-forge -c bioconda viroconstrictor
    ```

    You can change "mamba" for "conda" in the command above if you prefer to use conda or if mamba is not installed on your system.


## Download and installation from source

First start by cloning the repository and make sure you're on the latest released version of ViroConstrictor:

??? info inline end "Download a specific version"
    If you want to download a <u>**specific version**</u> of ViroConstrictor then please use the following command below where you replace `{VERSION}` with the version that you want.

    ```bash
    git clone -b {VERSION} --depth 1 https://github.com/RIVM-bioinformatics/ViroConstrictor.git; cd ViroConstrictor
    ```

```bash
# clone the repository from github and switch to the latest released version
git clone https://github.com/RIVM-bioinformatics/ViroConstrictor.git; cd ViroConstrictor; git checkout tags/$(git tag --sort=committerdate | tail -1) >> /dev/null
```


!!! tip "Make a new Conda environment before continuing"
    If you have Conda installed on your system, please create and activate a new environment before continuing.

    Use the following command to create and activate a new Conda environment named "ViroConstrictor" based on the environment-recipe we provide in the github-repository

    ```bash
    conda env create -f env.yml; conda activate ViroConstrictor
    ```
    The "ViroConstrictor" conda-environment should now be active

Once the ViroConstrictor repository is downloaded and you've created a conda environment to work in then you can install the downloaded package with pip:

```bash
# Install the downloaded package with pip
pip install .
```

ViroConstrictor should now be installed!
You can verify if installation was successful by typing `viroconstrictor --version` on the command-line, this should show the installed ViroConstrictor version.
