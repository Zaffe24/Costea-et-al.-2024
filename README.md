# Analyses from Costea et al. 2024

Repository containing the full set of scripts for the generation and visualization of the computational analyses in [Costea et al., 2024](https://doi.org/10.1101/2024.06.24.600391).

The scripts are organized according to the figures in the publication.

# Installation

To create a local conda environment with the required tools just use

`git clone https://github.com/Zaffe24/Costea-et-al.-2024`

`cd Costea-et-al.-2024`

`make all`

This usually takes a few minutes to create the conda environment.

# System requirements

We included a conda environment file, `environment.yml`, that lists all software dependencies.

`conda env create -f environment.yml`

For python packages, we also included a `requirements.txt` file.

`pip install -r requirements.txt`

# Further workflows

The single-cell alternative splicing analysis pipeline is part of a separate repository [AS_VASAseq_sc_pipeline](https://github.com/Zaffe24/AS_VASAseq_sc_pipeline).

# Sequencing data

The raw sequencing data can be accessed here: https://ega-archive.org/studies/EGAS50000000582




