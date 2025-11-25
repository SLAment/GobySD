# Characterizing the sex determination system of the two-spotted goby *Pomatoschistus flavescens*

Here you'll find the code associated with the genomic analyses of the manuscript:

Ament-Velásquez et al. (2025) "Female-biased sex ratios despite genetic sex determination in a marine fish with male-only parental care", biorxiv, []().

-----------

The pipelines were all desgined in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and depend on [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environments.

The standard with snakemake is that the pipeline is a file called `Snakefile`. However, I actually name them something else, like `mypipeline.smk` simply because it helps me to keep track of what pipeline is doing what other than based on it's path. Feel free to rename the files when you are working.

The directories in the repository are ordered to reflect the order of analyses in the paper. The folder `data` contains files used by multiple pipelines. The genome assemblies are too heavy for GitHub and will be provided elsewhere.

-----------

TODO:

- Remove relative paths in R scripts and in config files of some pipelines
- Describe folders here and mention what figures are produced where

----

Disclaimer: These scripts and files are provided "as is" and without any express or implied warranties, including, without limitation, the implied warranties of merchantability and fitness for a particular purpose.