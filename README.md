# Characterizing the sex determination system of the two-spotted goby *Pomatoschistus flavescens*

Here you'll find the code associated with the analyses of the manuscript:

Ament-Velásquez et al. (2025) "Female-biased sex ratios despite genetic sex determination across a climatic gradient in a marine fish", [biorxiv](https://www.biorxiv.org/content/10.64898/2026.02.16.706123v1.full).

## Data availability

- Smaller datasets and annotation files are provided in this repository. Larger files, such as the genome assemblies, GFF annotation files, and VCF variant files will be available in Zenodo (coming soon).
- The newly-generated genomic data for this study have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number [PRJEB103015](https://www.ebi.ac.uk/ena/browser/view/PRJEB103015) (coming soon).
- Previously produced ddRAD raw reads are available [here](http://www.ncbi.nlm.nih.gov/bioproject/1224570) (coming soon).

-----------

The bioinformatics pipelines were all designed in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and depend on [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environments.

The standard with snakemake is that the pipeline is a file called `Snakefile`. However, I actually name them something else, like `mypipeline.smk` simply because it helps me to keep track of what pipeline is doing what other than based on it's path. Feel free to rename the files when you are working, but instructions are provided with the existing names.

The order of the folders more or less reflects the order of figures in the paper. Most pipelines are designed to produce publication-ready figures, but some were further edited or assembled in [Inkscape](https://inkscape.org/).

The directories in the repository are ordered to reflect the order of analyses in the paper. The folder `data` contains files used by multiple pipelines.

-----------
TODO:

- Describe folders here and mention what figures are produced where
- Add link to the Zenodo repository

----

Disclaimer: These scripts and files are provided "as is" and without any express or implied warranties, including, without limitation, the implied warranties of merchantability and fitness for a particular purpose.