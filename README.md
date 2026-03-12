# Characterizing the sex determination system of the two-spotted goby *Pomatoschistus flavescens*

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

Here you'll find the code associated with the analyses of the manuscript:

Ament-Velásquez et al. (2025) "Female-biased sex ratios despite genetic sex determination across a climatic gradient in a marine fish", [biorxiv](https://www.biorxiv.org/content/10.64898/2026.02.16.706123v1.full).

## Data availability

- Smaller datasets and annotation files are provided in this repository. Larger files, such as the genome assemblies, GFF annotation files, and VCF variant files are available in [Zenodo](https://zenodo.org/records/18784321).
- The newly-generated genomic data for this study have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number [PRJEB103015](https://www.ebi.ac.uk/ena/browser/view/PRJEB103015) (coming soon).
- Previously produced ddRAD raw reads are available in National Center for Biotechnology Information (NCBI), [here](http://www.ncbi.nlm.nih.gov/bioproject/1224570) (coming soon).

**Looking for the genome assembly of our male fish (aka. the Norwegian individual)?**
Go to [Zenodo](https://zenodo.org/records/18784321)!

**Looking for the reference genome assembly of the two-spotted goby (aka. the French individual)?**
You can get it from NCBI (GenBank's accession number OZ251433.1) or from our [Zenodo](https://zenodo.org/records/18784321) repository.

**Looking for the genome annotations?**
Go to [Zenodo](https://zenodo.org/records/18784321)! But notice that some files are also in this repository, at `data/Annotation`.

**Looking for our repeat libraries?**
They are available here at `data/Annotation`.

-----------

The bioinformatics pipelines were all designed in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and depend on [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environments.

The standard with snakemake is that the pipeline is a file called `Snakefile`. However, I actually name them something else, like `mypipeline.smk` simply because it helps me to keep track of what pipeline is doing what other than based on its path. Feel free to rename the files when you are working, but instructions are provided with the existing names.

The order of the folders more or less reflects the order of figures in the paper. Most pipelines are designed to produce publication-ready figures, but some were further edited or assembled in [Inkscape](https://inkscape.org/).

The directories in the repository are ordered to reflect the order of analyses in the paper. The folder `data` contains files used by multiple pipelines.

----

Disclaimer: These scripts and files are provided "as is" and without any express or implied warranties, including, without limitation, the implied warranties of merchantability and fitness for a particular purpose.