
# README for Sex Determination of the two-spotted goby Pomatoschistus flavescens
--------------------------------------------------------------------------------

List of data files and scripts/code, organised by section of the project.

## 1. Analysis of field data of Adult Sex Ratio (ASR)
--------------------------------------------------
### Data: 
OSRdata.csv                                                                            # 2022 Field data of snorkeling observation of natural populations
sexrpheno.csv                                                                          # 2022 Field data of fish caught and phenotyped on land
Bergen_2019.txt, K_B_2007.txt, K_2010.txt, Tvarminne_2008.txt, beach_seines.txt        # Historical ASR data

### Code:
ASR_script.R      # Analyse 2022 ASR data and historical ASR data, produces corresponding figures and statistical tests. 


## 2. Identification of sex-specific loci using RADseq data, aligned to the reference genome.
-----------------------------------------------------------------------------------------

### Data:
ddRAD reads will be accessible at: http://www.ncbi.nlm.nih.gov/bioproject/1224570
reference genome will be accessible at: XXXX
popmap_hap1.txt      # Population map, indicates the sex and population of origin of each individual in the ddRAD dataset. For alignment to haplotype 1 of the reference genome.
popmap_hap2.txt      # Population map, indicates the sex and population of origin of each individual in the ddRAD dataset. For alignment to haplotype 2 of the reference genome.
popmap_radsex.tsv    # Population map used to identify the sex of individuals in the Radsex pipeline
TODO: link to vcf file repository

### Code:
index_flav.sh      # the reference genome is indexed before alignment of ddRAD reads
align_flav.sh      # Uses the ddRAD reads and the indexed reference genome and produces sorted alignment files, for haplotypes 1 and 2 of the reference genome.
gstacks_flav.sh    # Used the sorted alignment files and the population map and feeds them to the Stacks pipeline, first step gstacks to produce a catalog.
population_flav.sh # Uses the catalog and the population map to produce a VCF file.
assoc_flav.sh      # Filters the vcf, then uses plink to test for association of allele frequency and SNP missingness with sex
radsex_flav.sh     # Runs the Radsex pipeline on the ddRAD reads, finds sex specific markers, maps them to the reference genome
markers_list.R     # Extracts the list of markers associated with sex from the output of Stacks and Radsex. Also provides the calculation of the proportion of Heterozygous males and females for the sex markers in the ddRAD data.



## 3. Identification of the sex of eggs and larvae from genotype
--------------------------------------------------------------

### Data:
popmap_corrected           # Population map, indicates which population an nest the eggs and larvae come from, and in the case of adult Controls which sex they are.
general_stats_table.tsv    # Info on the percentage mapping to goby Reference genome of each sample
filesize.txt               # Info on number of reads of each sample
Sex_loci_RAD.txt           # List of sex associated markers identified by the RADsex pipeline
sexcalls2                  > TODO: link to vcf file repository

### Code:

sex_ID.R  # Identifies the genotypic sex of eggs and larvae, produces corresponding figures and Statistical tests.


