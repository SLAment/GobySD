#!/usr/bin/env Rscript

### GobyVCFmaf
#############################################################################
# Exploring the major allele distribution and genetic sex of gobies
# The data used here corresponds to the mapping of the ddRAD markers to the 
# P. minutus genome and predates our TH1 assembly and the reference fGobFla1.
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/10/17
# ++++++++++++++++++++++++++++++++++++++++++++++
# ============================
# Load the necessary libraries
# ============================
library(dplyr)
library(tidyr)
library(ggplot2)
library(vcfR) # The package to read vcf files!
library(stringr) # for str_remove
library(reshape2) # To fix dataframe with melt
library(ggtext) # to color labels
# ============================


# ============================
# Reading the data
# ============================
# IMPORTANT: set path of the repository
setwd("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Manuscripts/12_Goby/Repositories/GobySD/9_Others")

## File names
vcffile <- "data/sexcalls2.vcf.gz"
sexmarkersf <- "data/Sex_loci_RAD.txt"
sampleinfof <- "data/sampleinfo.txt"
popmapf <- "data/popmap_corrected.txt"
sexingmarkers <- read.table("../4_PlotMaleKmers/data/hap1_SexMarkers_selected_corrected.gff3", header = FALSE )

#### Read the files
vcf <- read.vcfR(vcffile, verbose = FALSE)
sexmarkers <- read.table(sexmarkersf, header = TRUE) # %>% mutate(Site = paste0(CHR, "_", "POS", "_"))
sampleinfo <- read.table(sampleinfof, header = TRUE) %>% mutate(sample_type = type) %>% select(-type) # Rename the type because it conflicts with the popmap later
popmap <- read.table(popmapf, header = TRUE) 

## Outputs
sexmarkersPCA <- "results/Sex_loci_PCA.png"
suppfigure <- "results/FigS1_MAF_Sex_loci_PCA.png"

# ============================
# Process vcf
# ============================

## Filter for the sex markers
# Extract the chromosome and position data from the VCF file
vcf_chrom <- getCHROM(vcf) # Chromosome information
vcf_pos <- getPOS(vcf)     # Position information

# Merge chromosome and position information into a data frame
vcf_sites <- data.frame(CHR = vcf_chrom, POS = vcf_pos)

# Create a logical vector to match rows in vcf_sites with those in sexmarkers
# This will be TRUE for rows in vcf that match chromosome (scaffold) and position in sexmarkers
matching_sites <- with(vcf_sites, (CHR %in% sexmarkers$CHR) & (POS %in% sexmarkers$POS))

# Now subset the VCF file to include only the matching sites
vcfsex <- vcf[matching_sites, ]
# ============================
# Get the final location of the sexy markers (TH1 coordinates)
# ============================
## Process the GFF3 into a more palatable table
gff2genenames <- function(gff){
  # version 3
  names(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  genenames <- c()
  IDs <- c()
  feattypes <- c()
  TRClist <- c()
  for(i in 1:nrow(gff)){
    attrib <- gff[i,9] %>% as.character() %>% strsplit(.,";") %>% .[[1]]
    name <- attrib[pmatch("Name=", attrib)] %>% strsplit(.,"=") %>% .[[1]] %>% .[2]
    genenames <- c(genenames, name)
    
    ID <- attrib[pmatch("ID=", attrib)] %>% strsplit(.,"=") %>% .[[1]] %>% .[2]
    IDs <- c(IDs, ID)
    
    TRC <- NA
    # feature_type <- gff[i, 3]
    if (grepl('TRC_', name)) {
      TRCstring <- strsplit(name,"_")[[1]]
      TRC <- paste0(TRCstring[1], '_', TRCstring[2])
      # feature_type <- 'satDNA'
    }
    TRClist <- c(TRClist, TRC)
    
    # if (grepl('simple_repeat', ID)) {
    #   feature_type <- 'simple_repeat'
    # } 
    # feattypes <- c(feattypes, feature_type)
  }
  gff <-  gff %>% select(seqid, start, end) %>% cbind(., "Name" = genenames, "ID" = IDs, "TRC" = TRClist) # , "type" = feattypes
  return(gff)
}

# Markers used for sexing
sexymarkers <- gff2genenames(sexingmarkers)

# ============================
# Calculate the MAF
# ============================
## This is a bcftools vcf
## Get the allele depth
ad <- extract.gt(vcfsex, element='AD') %>% data.frame()
ad <- tibble::rownames_to_column(ad, "Site") %>% data.frame()

# Reshape into long format
ad <- reshape2::melt(ad, id = c("Site"))
names(ad) <- c("Site", "Sample", "AD")
# Get the different allele depths as separate columns
ad <- ad %>% tidyr::separate(AD, c("Allele1", "Allele2", "Allele3"), sep = ",", remove = TRUE)
ad$Allele1 <- as.numeric(ad$Allele1)
ad$Allele2 <- as.numeric(ad$Allele2)
ad$Allele3 <- as.numeric(ad$Allele3)
# Often there is no third allele, so fix that to be 0 instead of NA
ad <- ad %>% mutate(Allele3 = case_when(is.na(Allele3) ~ 0, TRUE ~ Allele3)) # Fix the third allele

# Calculate their allele frequencies
ad <- ad %>% rowwise() %>% mutate(DP = Allele1 + Allele2 + Allele3) %>% mutate(freq1 = Allele1/DP, freq2 = Allele2/DP, freq3 = Allele3/DP) 

# Calculate the major (big) allele frequency 
MAFs <- ad %>% rowwise() %>% mutate(baf = max(Allele1, Allele2, Allele3)/DP) %>% data.frame() # rowwise makes the trick to do it per row
# # Add the genotypes (I didn't use this in the end)
# MAFs <- merge(MAFs, gt, by = c("Site", "Sample"))

## Add metadata
IDs <- str_remove(MAFs$Sample, "SORT.") %>% str_remove(., ".bam")
MAFs$ID <- as.numeric(IDs)
MAFs <- merge(MAFs, sampleinfo, by = "ID")
MAFs <- merge(MAFs, popmap, by = "ID")

## ==================
## Finding patterns
## ==================
## What samples have no coverage for any sex-linked marker?
# Make a new data frame for each sample to decide their genetic sex
SampleDF <- MAFs %>% 
  filter(!is.na(baf)) %>%  # remove missing data
  group_by(Sample) %>% summarise(totalDP = sum(DP), 
                                 avgBaf = mean(baf), # average major allele freq
                                 n = n()) %>% # count how many markers made it to the average
  mutate(genetic_sex = case_when(avgBaf > 0.9 ~ "Female", TRUE ~ "Male"))

# # Histogram of the total allele depth
# ggplot(SampleDF, aes(x = totalDP)) + geom_histogram()
# 
# # Global average major allele freq distribution
# ggplot(SampleDF, aes(x = avgBaf)) + geom_histogram() +
#   theme_bw() + xlab("Major allele frequency weighted by number of markers")
# # There are two very clear groups

## ==================
## Plot of the major allele frequency of the sexy markers
## ==================
# How does the major allele frequency look like if I divide them with genetic_sex
MAFs <- merge(MAFs, select(SampleDF, Sample, genetic_sex), by = "Sample")
MAFs <- MAFs %>% tidyr::separate(Site, c("Name", "pos_minutus"), remove = FALSE)

## Palette
femalecol <- "#FFD247"
malecol <- "#A482C4"

# Some markers were discarded
bad_markers <- c("scaffold2702_47682", "scaffold4193_6290", "scaffold709_69411", "scaffold709_69413")  # markers to highlight

label_fun <- function(value) {
  ifelse(value %in% bad_markers,
         paste0("<span style='color:red;'>", value, "</span>"),
         value) 
  }


### Major allele distribution of the markers for all samples
(majorallplot <- ggplot(MAFs, aes(baf, fill = genetic_sex)) + 
    geom_histogram() +
    facet_grid(genetic_sex ~ Site,
               labeller = labeller(Site = label_fun)) +
    xlab("Major allele frequency") +
    theme_bw() + 
    # scale_x_continuous(limits = c(0.5, 1), breaks = c(0.5, 0.75, 1)) +
    theme(legend.position = "none", 
          strip.background = element_blank(),
          axis.title = element_text(size = 14),
          strip.text.x = element_markdown(size = 9),  # key change to make labels red
          panel.grid.major = element_blank()) +
    scale_fill_manual(values= c("Female" = femalecol, "Male" = malecol))
)
# scaffold2702_47682 has a lot of missing data
# scaffold7908_69411 and scaffold709_69413 look strange, why is there a bump around 0.75? a paralog?
# scaffold4193_6290 is not exactly pretty but it goes more in the right direction

# Make a new data frame for each sample to decide their genetic sex
badSNPs <- c("scaffold2702_47682", "scaffold709_69411", "scaffold709_69413", "scaffold4193_6290")

## Controls
ggplot(MAFs %>% filter(type %in% c("male_control", "female_control")), aes(baf, fill = genetic_sex)) + 
  geom_histogram() +
  facet_grid(type~Site) + 
  xlab("Major allele frequency") +
  theme_bw() + theme(legend.position = "none") +
  scale_fill_manual(values= c("Female" = femalecol, "Male" = malecol))

MAFs %>% filter(type %in% c("male_control", "female_control")) %>% select(genetic_sex, type, ID) %>% distinct()
# Sample 745 should be male but looks genetically female

ggplot(MAFs %>% filter(Sample == "SORT.745.bam"), aes(baf)) + 
  geom_histogram() + facet_grid(Site~.) + xlim(0.5, 1.05)
# Mmm it looks indeed veeeery female.

## ==================
#### PCA of the samples with just the nice markers 
## ==================
minavgBaf <- 0.65
maxavgBaf <- 0.85

# Re-calculate avgBaf (average allele frequency of the major (big) allele)
SampleDFpretty <- MAFs %>% 
  filter(!is.na(baf)) %>%
  filter(!Site %in% badSNPs) %>% 
  group_by(Sample) %>% summarise(totalDP = sum(DP), 
                                 avgBaf = mean(baf), # average major allele freq
                                 n = n()) %>% # count how many markers made it to the average
  mutate(genetic_sex = case_when(avgBaf > maxavgBaf ~ "Female", 
                                 avgBaf < minavgBaf ~ "Male", 
                                 TRUE ~ "Undetermined"))

# SampleDFpretty %>% filter(Sample == "SORT.147.bam")
SampleDFpretty %>% filter(n < 5) # There are 3 samples with a single surviving marker... four samples have very low coverage

# Make a matrix of the allele frequencies
allele_freqs <- MAFs %>% select(Site, Sample, baf) %>%
  filter(!Site %in% badSNPs) %>% # The samples with missing markers are removed
  group_by(Site, Sample) %>%
  filter(!is.na(baf)) %>% # Remove rows with missing data
  tidyr::pivot_wider(names_from = Site, values_from = baf) # Spread the Sites into columns

# Check and remove any non-numeric columns just in case
allele_data <- allele_freqs[, sapply(allele_freqs, is.numeric)]
allele_data$Sample <- allele_freqs$Sample # put sample name back

# Remove samples and sites with missing data
allele_data_clean <- allele_data[complete.cases(allele_data), ]

# Perform PCA
pca_result <- prcomp(allele_data_clean %>% select(!Sample), center = TRUE, scale. = TRUE) 
setdiff(MAFs %>% select(Site, Sample, baf) %>% .$Sample, allele_data_clean$Sample) %>% length() # 8 samples were removed

# Get proportion of variance explained by each component
var_explained <- (pca_result$sdev)^2 / sum((pca_result$sdev)^2)
summary(pca_result)

# Extract the PCA scores for the first two components
pca_data <- as.data.frame(pca_result$x)
pca_data %>% dim() # 745 out of 768 samples

# Add sample IDs for labeling
pca_data$Sample <- allele_data_clean$Sample
pca_data$nums <- rownames(allele_data_clean)

# Add more information about the samples
pca_data <- merge(pca_data, SampleDFpretty, by = "Sample")
pca_data <- merge(pca_data, MAFs %>% select(Sample, type) %>% distinct(), by = "Sample")
pca_data <- pca_data %>% mutate(genetic_sex = case_when(avgBaf > maxavgBaf ~ "Female", 
                                                        avgBaf < minavgBaf ~ "Male", 
                                                        TRUE ~ "Undetermined"))

pca_data <- pca_data %>% mutate(Status = if_else(type %in% c("female_control", "male_control"), "Control", "Unknown"))

# The male that looks female
pca_data_outlier <- pca_data %>% filter(Sample == "SORT.745.bam")

## Plot the PCA 
(PCAp <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = totalDP, shape = genetic_sex)) +
    # ---- Some mixed samples were removed ----
  geom_vline(xintercept = -1.2, linetype = "dashed", color = "gray", linewidth = 1) +
  geom_vline(xintercept = 1.2, linetype = "dashed", color = "gray", linewidth = 1) +
    
    # ---- The samples ----
    theme_bw() + geom_point(alpha = 0.7, size = 3) +
    scale_colour_continuous(type = "viridis") +
    ggnewscale::new_scale("color") + # Geoms added to a plot after this function will use a new scale definition.
    geom_point(data = pca_data %>% filter(Status != "Unknown"), aes(colour = Status), size = 2) +
    geom_point(data = pca_data_outlier, aes(), colour = "black", shape = 10, size = 3) + # Outlier male control that looks like a female

    
    # ---- Add an annotation arrow and label ----
  annotate("segment",
           x = unique(pca_data_outlier$PC1) - 0.1,
           y = unique(pca_data_outlier$PC2) - 0.20,
           xend = unique(pca_data_outlier$PC1) - 0.02, 
           yend = unique(pca_data_outlier$PC2) - 0.05,
           arrow = arrow(length = unit(0.25, "cm")),
           colour = "black") +
    annotate("text",
             x = unique(pca_data_outlier$PC1) - 0.27,
             y = unique(pca_data_outlier$PC2) - 0.32,
             label = "Outlier\ncontrol",
             hjust = 0, size = 3.5) +
  # Labels
  xlab( paste0("PC1 (", round(var_explained[1],digits = 4), "% of variance)") ) +
  ylab( paste0("PC1 (", round(var_explained[2],digits = 4), "% of variance)") )  
)

# ============================
# Put them together for the final figure
# ============================

supfig <- cowplot::plot_grid(majorallplot, PCAp, 
                             nrow = 2, 
                             labels = c("a", "b"),
                             label_size = 20,
                             align = "h", axis = c("tblr"))

ggsave(suppfigure, plot = supfig, width = 12, height = 8)

