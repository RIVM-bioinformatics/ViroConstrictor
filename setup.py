import sys


if sys.version_info.major != 3 or sys.version_info.minor < 7:
    print("Error: you must execute setup.py using Python 3.7 or later")
    sys.exit(1)
    
from setuptools import setup, find_packages

exec(open("VAmp/version.py").read())

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="VAmp",
    description="Analysis of Viral Amplicon NGS data",
    author='Florian Zwagemaker, Dennis Schmitz',
    author_email='rivm-bioinformatics@rivm.nl',
    license='AGPLv3',
    version=__version__,
    packages=find_packages(),
    scripts=[
        'VAmp/workflow/workflow.smk',
        'VAmp/workflow/directories.py'],
    package_data={'VAmp': ['workflow/envs/*', 'workflow/scripts/*', 'workflow/files/*']},
    install_requires=[
        'gnureadline>=8.0.0',
        'biopython>=1.78',
        'snakemake>=6.0.5'
    ],
    entry_points={"console_scripts": [
        'vamp = VAmp.VAmp:main',
        'VAmp = VAmp.VAmp:main']},
    keywords=[],
    include_package_data=True,
    zip_safe=False
)