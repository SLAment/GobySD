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
library(tidyr)
# ============================

# ============================
# Input
# ============================
# Input
covdf <- read.table(snakemake@input$covwin, header = TRUE)
gffsatDNAraw <- read.table(snakemake@input$gff, header = FALSE)
gffamhr2y <- read.table(snakemake@input$gene) %>% filter(V3 == "exon")
STARTSD <- snakemake@params$start

# Output
SDfile <- snakemake@output$chrs

# ============================
# Processing
# ============================

covdf <- covdf %>% mutate(Midpos = ((end - start)/2) + start)

# Calculate an average coverage per sample but using the part of the scaffold 
# before the SD region
meancovs <- covdf %>% filter(end < STARTSD) %>% 
  group_by(Sample) %>% summarise(mean = mean(Coverage), sdcov = sd(Coverage) )

# Add the geographic population
covdf <- covdf %>%
  mutate(Population = case_when(
    stringr::str_starts(Sample, "PflaHEL") ~ "Helligvaer",
    stringr::str_starts(Sample, "PflaKGB") ~ "Kristineberg"
  ))

# Calculate a normalized coverage 
covdf_norm <- covdf %>% 
  left_join(meancovs, by = "Sample") %>%
  mutate(NormCoverage = Coverage / mean)

# covdf[apply(covdf, 1, function(x) any(is.na(x))), ] # No window has 0 coverage

# A plot comparing the difference in coverage
covdf_norm_diff <- covdf_norm %>% 
  group_by(seqid, Midpos) %>%
  summarise(
    male_cov = mean(NormCoverage[Sex == "Male"], na.rm = TRUE),
    female_cov = mean(NormCoverage[Sex == "Female"], na.rm = TRUE)
  ) %>%
  mutate(diff = female_cov - male_cov)

# Similar but emphasis on the coverage itself, not the difference
cov_summary <- covdf_norm %>%
  group_by(seqid, start, end, Midpos, Sex) %>%
  summarise(mean_cov = mean(NormCoverage, na.rm = TRUE)) %>%
  pivot_wider(names_from = Sex, values_from = mean_cov, names_prefix = "cov_")

filtered_windows <- cov_summary %>%
  filter(cov_Female < 0.15,
         cov_Male >= 0.3,
         cov_Male <= 1.3)

# Prepare annotation of genes
names(gffamhr2y) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# names(gffgenes) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

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

gffsatDNA <- gff2satDNA(gffsatDNAraw) %>% filter(seqid == "h1tg000060l")

# Define custom color mapping
my_colors <- c(
  "TRC_1"   = "#FE8776",
  "TRC_3"   = "#3DBEAC",
  "sf10__TRC_75"  = "#AE123A",
  "sf11__TRC_42"  = "#FF00FF",
  "TRC_105" = "#F8AF50" )

# ============================
# Plotting
# ============================
## Palette
femalecol <- "#FFD247"
malecol <- "#A482C4"

# ## Distribution of coverage
# (covdensplot <- ggplot(covdf %>%
#                          filter(Coverage < 80, Coverage > 20),
#                        aes(x = Coverage, fill = Sex)) +
#     geom_density() +
#     facet_grid(seqid ~ Sample) +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank()) +
#     scale_fill_manual(values= c("Female" = femalecol, "Male" = malecol))
# )
# 
# ggsave(filename = densityfile,
#        plot = covdensplot,
#        width = 7, height = 7)

( SDhap1 <- ggplot(covdf_norm, 
                   aes(x = Midpos, y = NormCoverage, line = Population, colour = Sex)) +
    geom_hline(yintercept=1, colour = "gray30") +
    # geom_point(data = filtered_windows %>% mutate(Population = "Helligvaer"),
    #            aes(y = cov_Female),
    #            color = "red", size = 1, alpha = 0.7) +
    geom_line(alpha = 0.8, linewidth = 0.6, aes(linetype = Population)) +
    theme_bw() +
    scale_colour_manual(values= c("Female" = femalecol, "Male" = malecol)) +
    xlab("Coordinate of the Y haplotype (bp)") +
    ylab("Normalized depth") +
    # coord_cartesian(ylim = c(-0.5, 5), x = c(STARTSD, 4e+06)) +
    coord_cartesian(ylim = c(-0.3, 3), x = c(STARTSD, 4e+06)) +
    
    # ## genes
    # geom_rect(data = gffgenes,
    #           aes(xmin = start, xmax = end, ymin = -0.6, ymax = -0.5),
    #           colour = "black",
    #           inherit.aes = FALSE,
    #           alpha = 0.1,
    #           show.legend = FALSE) +
    
    ## amhr2y
    geom_rect(data = gffamhr2y,
              aes(xmin = start, xmax = end, ymin = -0.3, ymax = -0.2),
              colour = "black",
              inherit.aes = FALSE,
              alpha = 0.5,
              show.legend = FALSE) +
    
    ## Add the satDNA annotation
    geom_rect(data = gffsatDNA,
              aes(xmin = start, xmax = end, ymin = -0.3, ymax = -0.2,
                  fill = Name),
              inherit.aes = FALSE,
              alpha = 0.5,
              show.legend = FALSE) +
    scale_fill_manual(
      values = my_colors,
      na.value = "grey40")   ) # all others not in my_colors will be grey

ggsave(filename = SDfile, 
       plot = SDhap1,
       width = 10, height = 5)




## Not used
# ============================
# Difference in coverage
# ============================
ggplot(covdf_norm_diff, 
       aes(x = Midpos, y = diff)) +
  geom_hline(yintercept=0.5, colour = "gray") +
  geom_hline(yintercept=-0.5, colour = "gray") +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  theme_bw() +
  xlab("Chromosome coordinate (bp)") +
  ylab("Coverage difference between males and females") +
  # geom_point(data = filtered_windows,
  #            aes(y = cov_Female),
  #            color = "red", size = 1, alpha = 0.7) +
  
  ## Add the satDNA annotation
  geom_rect(data = gffsatDNA,
            aes(xmin = start, xmax = end, ymin = 0.1, ymax = 0,
                fill = Name),
            inherit.aes = FALSE,
            alpha = 0.5,
            show.legend = FALSE) +
  ## amhr2y
  geom_rect(data = gffamhr2y,
            aes(xmin = start, xmax = end, ymin = 0.1, ymax = 0),
            colour = "red",
            inherit.aes = FALSE,
            alpha = 0.5,
            show.legend = FALSE) +
  scale_fill_manual(
    values = my_colors,
    na.value = "grey40") 


