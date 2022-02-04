# Installation

ViroConstrictor is only available on Linux (or Linux-based) operating systems. MacOS may also work but is not tested.
ViroConstrictor will *not work* on Windows.  

***Download and installation instructions through Conda will be made available as soon as possible.***  
Until then, please perform manual installation as described on this page below.

## Download

Use the following command to download the <u>**latest**</u> release of ViroConstrictor and move to the newly downloaded `ViroConstrictor/` directory:

```bash
git clone https://github.com/RIVM-bioinformatics/ViroConstrictor.git; cd ViroConstrictor; git checkout tags/$(git tag --sort=committerdate | tail -1) >> /dev/null
```

!!! tip  
    If you want to download a <u>**specific version**</u> of ViroConstrictor then please use the following command below where you replace `{VERSION}` with the version that you want.

    ```bash
    git clone --depth 1 -b {VERSION} https://github.com/RIVM-bioinformatics/ViroConstrictor.git; cd ViroConstrictor
    ```

## Installation

**Before you install ViroConstrictor, make sure [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) is installed on your system and functioning properly!**

1. Create the required conda-environment and install the necessary dependencies.  
    Copy and paste the code-snippet below to create the new conda-environment and directly activate it.  
    ```bash
    conda create --name ViroConstrictor -c conda-forge mamba python=3.7 -y; 
    conda activate ViroConstrictor; mamba env update -f mamba-env.yaml
    ```

    ??? info "Click here if that code-snippet did not work for you"
        if you get errors from Conda, you can also try the following command(s):
        ```bash
        conda env create -f env.yml && conda activate ViroConstrictor
        ```
    The "ViroConstrictor" environment should now be active  

2. You can now actually install ViroConstrictor on your system with:  
    ```
    pip install .
    ```

    ??? info "Click here if that command did not work for you"
        If you get errors from `pip`, you can also try the following command: 
        ```bash
        python setup.py install
        ```

**ViroConstrictor is now installed!**  
You can start ViroConstrictor from anywhere on your system as long as the ViroConstrictor conda-environment is active.  
You can also use the ViroConstrictor pipeline in a different conda-environment as long as the software dependencies match.

You can verify if ViroConstrictor was properly installed by typing `viroconstrictor -v`, which should return the installed version. Or by typing `viroconstrictor -h` which should return the help document.