# Analyses related to adult sex ratio in populations of the two-spotted goby Pomatoschistus flavescens
--------------------------------------------------------------------------------

List of data files and scripts/code, organised by section of the project.

## 1. Analysis of field data of Adult Sex Ratio (ASR)
--------------------------------------------------

### List of datasets: 
<ins>**OSRdata.csv:**</ins> 2022 Field data of snorkeling observation of natural populations

<ins>**sexrpheno.csv:**</ins> 2022 Field data of fish caught and phenotyped on land

<ins>**Bergen_2019.txt, K_B_2007.txt, K_2010.txt, Tvarminne_2008.txt, beach_seines.txt:**</ins> Historical ASR data

### Code:
<ins>**ASR_script.R:**</ins>  Analyse 2022 ASR data and historical ASR data, produces corresponding figures and statistical tests. 


## 2. Identification of sex-specific loci using RADseq data, aligned to the reference genome.
-----------------------------------------------------------------------------------------

### External data:
ddRAD reads will be accessible at: http://www.ncbi.nlm.nih.gov/bioproject/1224570 

reference genome will be accessible at: XX Zenodo repository to be updated XX

vcf files will be accessible at: XX Dryad repository to be updated XX

### List of datasets: 
<ins>**popmap_hap1.txt :**</ins>     Population map, indicates the sex and population of origin of each individual in the ddRAD dataset. For alignment to haplotype 1 of the reference genome.

<ins>**popmap_hap2.txt:**</ins>      Population map, indicates the sex and population of origin of each individual in the ddRAD dataset. For alignment to haplotype 2 of the reference genome.

<ins>**popmap_radsex.tsv:**</ins>    Population map used to identify the sex of individuals in the Radsex pipeline

### Code:
<ins>**index_flav.sh:**</ins>      The reference genome is indexed before alignment of ddRAD reads.

<ins>**align_flav.sh:**</ins>      Uses the ddRAD reads and the indexed reference genome and produces sorted alignment files, for haplotypes 1 and 2 of the reference genome.

<ins>**gstacks_flav.sh:**</ins>   Used the sorted alignment files and the population map and feeds them to the Stacks pipeline, first step gstacks to produce a catalog.

<ins>**population_flav.sh:**</ins> Uses the catalog and the population map to produce a VCF file.

<ins>**assoc_flav.sh:**</ins>     Filters the vcf, then uses plink to test for association of allele frequency and SNP missingness with sex.

<ins>**radsex_flav.sh:**</ins>     Runs the Radsex pipeline on the ddRAD reads, finds sex specific markers, maps them to the reference genome.

<ins>**markers_list.R:**</ins>     Extracts the list of markers associated with sex from the output of Stacks and Radsex. Also provides the calculation of the proportion of Heterozygous males and females for the sex markers in the ddRAD data.



## 3. Identification of the sex of eggs and larvae from genotype
--------------------------------------------------------------

### List of datasets:
<ins>**popmap_corrected:**</ins>           Population map, indicates which population an nest the eggs and larvae come from, and in the case of adult Controls which sex they are.

<ins>**general_stats_table.tsv:**</ins>    Info on the percentage mapping to goby Reference genome of each sample.

<ins>**filesize.txt:**</ins>               Info on number of reads of each sample.

<ins>**Sex_loci_RAD.txt:**</ins>           List of sex associated markers identified by the RADsex pipeline.

<ins>**sexcalls2.vcf.gz :**</ins>          VCF file of genotyped eggs and larvae aligned to the reference genome, in folder 09_Others/data

### Code:

<ins>**sex_ID.R:**</ins>  Identifies the genotypic sex of eggs and larvae, produces corresponding figures and Statistical tests.




