#!/usr/bin/env Rscript

### GobyIlluCovDist: Coverage of the Illumina goby data mapped to the French reference
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/10/08
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(cowplot)
library(stringr)
# ============================

# ============================
# Input
# ============================
# Input
covdf <- read.table(snakemake@input$covwin, header = TRUE)
gffsatDNAraw <- read.table(snakemake@input$gff, header = FALSE)

# Output
densityfile <- snakemake@output$density
chrsfile <- snakemake@output$chrs

# ## Coverage BED files
# covdf <- read.table("/Users/lorena/Dropbox/VRwork/Analyses/7_Goby/Dardel/09_IlluCoverageOnFrench/results/Goby_covwin30kb.bed", header = TRUE)
# gffsatDNAraw <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/IFB/07_RepeatMaskerAgain/results/fGobFla1-GobyTideCluster_v1.00.gff3", header = FALSE)
# SDchr <- "OZ251425.1"
# 
# # Output
# densityfile <- "/Users/lorena/Dropbox/VRwork/Analyses/7_Goby/Dardel/09_IlluCoverageOnFrench/results/Coverage_density_chrs.png"
# chrsfile <- "/Users/lorena/Dropbox/VRwork/Analyses/7_Goby/Dardel/09_IlluCoverageOnFrench/results/Coverage_scan_chrs.png"

# ============================
# Processing
# ============================

covdf <- covdf %>% mutate(Midpos = ((end - start)/2) + start)

# Change names to match the chromosomes
covdf <- covdf %>%
  mutate(seqid = recode(seqid,
                        "OZ251410.1" = "Chr1",
                        "OZ251411.1" = "Chr2",
                        "OZ251412.1" = "Chr3",
                        "OZ251415.1" = "Chr6",
                        "OZ251413.1" = "Chr4",
                        "OZ251414.1" = "Chr5",
                        "OZ251416.1" = "Chr7",
                        "OZ251417.1" = "Chr8",
                        "OZ251418.1" = "Chr9",
                        "OZ251419.1" = "Chr10",
                        "OZ251420.1" = "Chr11",
                        "OZ251421.1" = "Chr12",
                        "OZ251422.1" = "Chr13",
                        "OZ251423.1" = "Chr14",
                        "OZ251424.1" = "Chr15",
                        "OZ251425.1" = "Chr16",
                        "OZ251426.1" = "Chr17",
                        "OZ251427.1" = "Chr18",
                        "OZ251428.1" = "Chr19",
                        "OZ251429.1" = "Chr20",
                        "OZ251430.1" = "Chr21",
                        "OZ251431.1" = "Chr22",
                        "OZ251432.1" = "Chr23",
                        "OZ251433.1" = "Mt"))

# Sort the chromosomes in a natural order
covdf <- covdf %>%
  mutate(seqid = factor(seqid,
                        levels = c(paste0("Chr", 1:23), "Mt") ))  # or 1:13, or 1:30, etc.


# Calculate an average coverage
meancovs <- covdf %>% group_by(Sample) %>% summarise(mean = mean(Coverage), sdcov = sd(Coverage) )

# Add the geographic population
covdf <- covdf %>%
  mutate(Population = case_when(
    stringr::str_starts(Sample, "PflaHEL") ~ "Helligvaer",
    stringr::str_starts(Sample, "PflaKGB") ~ "Kristineberg"
  ))


# Calculate a normalized coverage per sample
covdf_norm <- covdf %>%
  left_join(meancovs, by = "Sample") %>%
  mutate(NormCoverage = Coverage / mean)

# ============================
# Prepare annotation
# ============================
## Process the GFF3 into a more palatable table, vectorized version: satDNA
gff2satDNA <- function(gff){
  names(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Extract Name and ID using regex
  gff <- gff %>%
    mutate(
      Name = str_extract(attributes, "Name=[^;]+") %>% str_replace("Name=", ""),
      ID = str_extract(attributes, "ID=[^;]+") %>% str_replace("ID=", ""),
      TRC = ifelse(str_detect(Name, "TRC_"),
                   str_match(Name, "(TRC_\\d+)")[,2],
                   NA)) %>% 
    select(seqid, start, end, Name, ID, TRC)
  
  return(gff)
}

gffsatDNA <- gff2satDNA(gffsatDNAraw) %>% 
  mutate(seqid = recode(seqid,
                        "OZ251410.1" = "Chr1",
                        "OZ251411.1" = "Chr2",
                        "OZ251412.1" = "Chr3",
                        "OZ251415.1" = "Chr6",
                        "OZ251413.1" = "Chr4",
                        "OZ251414.1" = "Chr5",
                        "OZ251416.1" = "Chr7",
                        "OZ251417.1" = "Chr8",
                        "OZ251418.1" = "Chr9",
                        "OZ251419.1" = "Chr10",
                        "OZ251420.1" = "Chr11",
                        "OZ251421.1" = "Chr12",
                        "OZ251422.1" = "Chr13",
                        "OZ251423.1" = "Chr14",
                        "OZ251424.1" = "Chr15",
                        "OZ251425.1" = "Chr16",
                        "OZ251426.1" = "Chr17",
                        "OZ251427.1" = "Chr18",
                        "OZ251428.1" = "Chr19",
                        "OZ251429.1" = "Chr20",
                        "OZ251430.1" = "Chr21",
                        "OZ251431.1" = "Chr22",
                        "OZ251432.1" = "Chr23",
                        "OZ251433.1" = "Mt"))

# Define custom color mapping
my_colors <- c(
  "TRC_1"   = "#FE8776",
  "TRC_3"   = "#3DBEAC",
  "sf10__TRC_75"  = "slateblue3",
  "sf11__TRC_42"  = "#FF00FF",
  "TRC_105" = "#AE123A" )

# ============================
# Plotting
# ============================
## Palette
femalecol <- "#FFD247"
malecol <- "#A482C4"

## Distribution of coverage
(covdensplot <- ggplot(covdf %>% 
                         filter(Coverage < 80, Coverage > 20), 
                       aes(x = Coverage, fill = Sex)) + 
    geom_density() +
    facet_grid(seqid ~ Sample) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(values= c("Female" = femalecol, "Male" = malecol))
)

ggsave(filename = densityfile, 
       plot = covdensplot,
       width = 7, height = 14)

# Plot the sex chr coverage in particular
( fullchrs <- ggplot(covdf_norm %>% filter(seqid %in% c("Chr15", "Chr16", "Chr17")), 
       aes(x = Midpos, y = NormCoverage, line = Population, colour = Sex)) +
  geom_hline(yintercept=1, colour = "gray30") +
  geom_line(alpha = 0.8, linewidth = 0.6, aes(linetype = Population)) +
  facet_grid(seqid ~ .) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3)) +
  scale_colour_manual(values= c("Female" = femalecol, "Male" = malecol)) +
  xlab("Chromosome coordinate (bp)") +
  theme(axis.title.x = element_blank()) +
  ylab("Normalized depth") +
  ## Add the satDNA annotation
  geom_rect(data = gffsatDNA %>% filter(seqid %in% c("Chr15", "Chr16", "Chr17")),
            aes(xmin = start, xmax = end, ymin = 0.1, ymax = 0,
                fill = Name),
            inherit.aes = FALSE,
            alpha = 0.5,
            show.legend = FALSE) +
  scale_fill_manual(
    values = my_colors,
    na.value = "grey40")   # all others not in my_colors will be grey
)

SDchr16 <- ggplot(covdf_norm %>% filter(seqid %in% c("Chr16")), 
       aes(x = Midpos, y = NormCoverage, line = Population, colour = Sex)) +
  geom_hline(yintercept=1, colour = "gray30") +
  geom_line(alpha = 0.8, linewidth = 0.6, aes(linetype = Population)) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3), xlim = c(2.85e+07, 3.43e+07)) +
  facet_grid(seqid ~ .) +
  scale_colour_manual(values= c("Female" = femalecol, "Male" = malecol)) +
  xlab("Chromosome coordinate (bp)") +
  ylab("Normalized depth") +
  theme(legend.position="none") +
  ## Add the satDNA annotation
  geom_rect(data = gffsatDNA %>% filter(seqid %in% c("Chr16")),
            aes(xmin = start, xmax = end, ymin = 0.1, ymax = 0,
                fill = Name),
            inherit.aes = FALSE,
            alpha = 0.5,
            show.legend = FALSE) +
  scale_fill_manual(
    values = my_colors,
    na.value = "grey40")   # all others not in my_colors will be grey


chrcomparisons <- cowplot::plot_grid(fullchrs, SDchr16, nrow=2, 
                   labels=c('a', 'b'), 
                   rel_heights = c(3, 1),
                   align = "vh", axis = "lrtb") 

ggsave(filename = chrsfile, 
       plot = chrcomparisons,
       width = 7, height = 10)
