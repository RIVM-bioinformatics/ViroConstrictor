import sys

from packaging import version as vv
from setuptools import find_packages, setup

from ViroConstrictor import __version__

if sys.version_info.major != 3 or sys.version_info.minor < 10:
    print("Error: you must execute setup.py using Python 3.10 or later")
    sys.exit(1)

try:
    import conda
except SystemError:
    sys.exit(
        """
Error: conda could not be accessed.
Please make sure conda is installed and functioning properly before installing ViroConstrictor
"""
    )

try:
    import snakemake
except SystemError:
    sys.exit(
        """
Error: SnakeMake could not be accessed.
Please make sure SnakeMake is installed properly before installing ViroConstrictor
"""
    )

if vv.parse(snakemake.__version__) < vv.parse("7.15.2"):
    sys.exit(
        f"""
The installed SnakeMake version is older than the minimally required version:

Installed SnakeMake version: {snakemake.__version__}
Required SnakeMake version: 7.15.2 or later

Please update SnakeMake to a supported version and try again
"""
    )

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="ViroConstrictor",
    description="Analysis of Viral Amplicon NGS data",
    author="Florian Zwagemaker, Dennis Schmitz, Karim Hajji, Annelies kroneman",
    author_email="ids-bioinformatics@rivm.nl",
    license="AGPLv3",
    version=__version__,
    packages=find_packages(),
    python_requires=">=3.10",
    scripts=[
        "ViroConstrictor/workflow/workflow.smk",
        "ViroConstrictor/workflow/match_ref.smk",
        "ViroConstrictor/workflow/directories.py",
        "ViroConstrictor/workflow/presets.py",
    ],
    package_data={
        "ViroConstrictor": [
            "workflow/envs/*",
            "workflow/scripts/*",
            "workflow/files/*",
            "workflow/scripts/match_ref/*",
        ]
    },
    install_requires=[
        "urllib3",
        "biopython>=1.79",
        "drmaa==0.7.9",
        "fpdf2",
        "pandas>=1.4.2",
        "openpyxl",
        "pyyaml==6.0",
        "rich==13.*",
        "AminoExtract==0.3.1",
    ],
    entry_points={
        "console_scripts": [
            "ViroConstrictor = ViroConstrictor.__main__:main",
            "viroconstrictor = ViroConstrictor.__main__:main",
            "viroConstrictor = ViroConstrictor.__main__:main",
            "Viroconstrictor = ViroConstrictor.__main__:main",
        ]
    },
    keywords=[],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
