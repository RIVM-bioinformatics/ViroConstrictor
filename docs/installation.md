# Installation

ViroConstrictor is available only on Linux (or Linux-based) operating systems. MacOS may also work, but it has not been tested. ViroConstrictor will *not work* on Windows.  

## Prerequisites

ViroConstrictor requires the following software to be installed on your system:

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://mamba.readthedocs.io/en/latest/installation.html)

### Recommended (Optional)

- [Apptainer](https://apptainer.org/)

## Installation with Conda

ViroConstrictor releases are published on Bioconda to allow for easy installation.

If you already have [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) installed, please use the command below:
```bash
mamba install -c conda-forge -c bioconda viroconstrictor
```

??? info "Click here if you don't have mamba installed in your system already"
    ```bash
    conda install -c conda-forge -c bioconda viroconstrictor
    ```

Installation through Mamba/Conda is the recommended method for installing ViroConstrictor.

!!! tip "Use a Dedicated Conda Environment for ViroConstrictor"
    It is recommended to install and use ViroConstrictor in a dedicated Conda environment. This will help ensure that there are no conflicting dependencies.

    You can create a dedicated environment and install ViroConstrictor in it using a single command with either Mamba or Conda.

    The following command will use **Mamba** to create a new environment named `pipeline_env` with ViroConstrictor installed:
    ```bash
    mamba create --name pipeline_env -c conda-forge -c bioconda viroconstrictor
    ```

    You can replace "mamba" with "conda" in the command above if you prefer to use Conda or if Mamba is not installed on your system.



## Download and installation from source

First, clone the repository and make sure you're on the latest released version of ViroConstrictor:

??? info inline end "Download a specific version"
    If you want to download a <u>specific version</u> of ViroConstrictor, use the following command, replacing `{VERSION}` with your desired version.

    ```bash
    git clone -b {VERSION} --depth 1 https://github.com/RIVM-bioinformatics/ViroConstrictor.git; cd ViroConstrictor
    ```

```bash
# Clone the repository from GitHub and switch to the latest released version
git clone https://github.com/RIVM-bioinformatics/ViroConstrictor.git; cd ViroConstrictor; git checkout tags/$(git tag --sort=committerdate | tail -1) >> /dev/null
```


!!! tip "Create a New Conda Environment Before Continuing"
    Use the following command to create and activate a new Conda environment named `ViroConstrictor` based on the environment recipe provided in the GitHub repository:

    ```bash
    conda env create -f env.yml; conda activate ViroConstrictor
    ```
    The "ViroConstrictor" Conda environment should now be active

Once you have downloaded the ViroConstrictor repository and created a Conda environment to work in, install the package with pip:

```bash
# Install the downloaded package with pip
pip install .
```

ViroConstrictor should now be installed!
You can verify if the installation was successful by typing `viroconstrictor --version` on the command line; this should display the installed ViroConstrictor version.