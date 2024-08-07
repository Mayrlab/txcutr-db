# Overview
[![Snakemake-Validate](https://github.com/Mayrlab/txcutr-db/actions/workflows/snakemake-validate.yml/badge.svg)](https://github.com/Mayrlab/txcutr-db/actions/workflows/snakemake-validate.yml)
[![DOI](https://zenodo.org/badge/414776529.svg)](https://zenodo.org/badge/latestdoi/414776529)


This repository provides a Snakemake pipeline for generating the target files for use in 
[the scUTRquant pipeline](https://Mayrlab.github.io/scUTRquant). This is provided as a record
of how we generated truncated transcriptomes for [the scUTRquant manuscript](https://www.biorxiv.org/content/10.1101/2021.11.22.469635v1) 
and an example of how to use the Bioconductor package [`txcutr`](https://bioconductor.org/packages/release/bioc/html/txcutr.html). 

Please note that, while the pipeline does provide some flexibility, it was implemented with the limited
scope of `mm10` and `hg38` annotations from Ensembl and GENCODE. For example, it must be modified
in ordered to generate correct FASTA files for `mm39` or `hg38` references.

# Setup
## Prerequisites
- Snakemake >= 5.11
- Conda/Mamba
- (optional) CellRanger

This should be compatible with Linux and MacOS systems. If Conda is not already installed, we recommend 
installing [Miniforge](https://github.com/conda-forge/miniforge#miniforge).

## Installation

```bash
git clone https://github.com/Mayrlab/txcutr-db.git
```

## Configuration
Please edit the `config.yaml` file to provide a `tmpdir` specific to your system. If you wish 
to use the GTF filtering provided by CellRanger, also specify the path to CellRanger for your system.

## Usage
The `rule all:` in the Snakefile contains specifications for several variants that were used in the scUTRquant 
manuscript. One likely does not want to generate all of these. Instead, a single variant can be "requested" at
the commandline. Since the `kallisto` index (.kdx file) is the last output, that is what should be specified:

```bash
snakemake --use-conda homo_sapiens/gencode.v38.annotation.pc.txcutr.w500.kdx
```

This would use the GENCODE v38 annotation, filtered for only protein-coding transcripts (`.pc`) with validated
3' ends, and truncated to 500 nts (`.w500`). The default merge table (TSV) will use a 200 nt merge distance.

# Notes

The `txcutr` step is computationally demanding. For example, in an HPC setting, we have it configured to 
run with 20 cores and 4 GB/core, which takes about 30 mins. 

Be aware that some rules include `thread` and `resources` specifications that are used by Snakemake cluster 
profiles. Please adjust accordingly (e.g., not all cluster configurations interpret the `mem_mb` parameter 
as *per core*)!
