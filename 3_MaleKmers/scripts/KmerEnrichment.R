#!/usr/bin/env Rscript

### KmerEnrichment: Looking for SD candidates
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/10/06
# ++++++++++++++++++++++++++++++++++++++++++++++
# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(tidyverse)
library(ggplot2)
library(purrr) # for fancy enrichment
library(ggrepel)
library(cowplot)
# ============================

# ============================
# Input files
# ============================
# Input
kmercoverageh1 <- read.table(snakemake@input$hap1)
kmercoverageh2 <- read.table(snakemake@input$hap2)
# Output
plotfile <- snakemake@output$plot

# # Input
# kmercoverageh1 <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/Dardel/04_MaleKmers/reports/male_coverage_hap1.txt", header = FALSE)
# kmercoverageh2 <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/Dardel/04_MaleKmers/reports/male_coverage_hap2.txt", header = FALSE)
# # Output
# plotfile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/Dardel/04_MaleKmers/results/FigS4_EnrichmentKmer.png"


# ============================
# Prepare data frame
# ============================
kmercoverage <- rbind(kmercoverageh1 %>% mutate(Haplotype = "Haplotype 1"),
                      kmercoverageh2 %>% mutate(Haplotype = "Haplotype 2"))

names(kmercoverage) <- c("Scaffold", "Start", "End", "Coverage", "Haplotype")

# define a window as “male-associated” if its Coverage > 0.05
kmercoverage <- kmercoverage %>% mutate(Male_associated = Coverage > 0.05)

# ============================
# Male-associated k-mer coverage distribution
# ============================

kmerdist <- ggplot(kmercoverage, aes(x = Coverage)) + geom_histogram() +
  theme_classic() +
  facet_grid(Haplotype ~ .) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "red") 
  

# ============================
# Male-associated k-mer coverage enrichment
# ============================
# Count how many windows per scaffold exceed the threshold, and what proportion that represents:

scaffold_summary <- kmercoverage %>%
  group_by(Haplotype, Scaffold) %>%
  summarise(
    n_windows = n(),
    n_male = sum(Male_associated),
    prop_male = n_male / n_windows
  ) %>%
  arrange(desc(prop_male))

enrichment <- kmercoverage %>%
  # group_map(~{ ... }) is the same as group_map(function(.x, .y) { ... })
  # .x = the data for one group (a tibble containing all rows for that group)
  # .y = the group keys (e.g. Haplotype = "Haplotype 1")
  group_by(Haplotype) %>%
  group_map(~{. # That tilde introduces an anonymous function in tidyverse syntax
    hap_data <- .x
    scaffolds <- unique(hap_data$Scaffold)
    
    # combine the returned tibbles row-wise
    purrr::map_dfr(scaffolds, function(scaf) { 
      scaf_mask <- hap_data$Scaffold == scaf
      
      tab <- table(
        In_scaffold = scaf_mask,
        Male_associated = hap_data$Male_associated
      )
      
      # Some scaffolds may have 0/0 categories
      if (any(tab == 0)) {
        pval <- NA
      } else {
        pval <- fisher.test(tab)$p.value
      }
      
      tibble(
        Scaffold = scaf,
        p_value = pval
      )
    })
  }) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Add the haplotype back
enrichment <- enrichment %>%
  left_join(kmercoverage %>% select(Scaffold, Haplotype) %>% distinct(), by = "Scaffold")

# Make the final enrichment table with all info
final_enrichment <- kmercoverage %>%
  group_by(Haplotype, Scaffold) %>%
  summarise(
    n_windows = n(),
    n_male = sum(Male_associated),
    prop_male = n_male / n_windows,
    .groups = "drop"
  ) %>%
  left_join(enrichment, by = c("Haplotype", "Scaffold"))

sig_scaffolds <- final_enrichment %>% filter(p_adj < 0.05, prop_male > 0.05)

(enrichplot <- ggplot(final_enrichment, aes(x = prop_male, y = -log10(p_adj), color = Haplotype)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = sig_scaffolds,
    aes(label = Scaffold),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  theme_bw() +
  labs(
    x = "Proportion of male-associated windows",
    y = "-log10(FDR-adjusted p-value)") )


# ============================
# Put together
# ============================

(FinalFigure <- cowplot::plot_grid(kmerdist, 
                                   enrichplot, 
                                   ncol = 2, 
                                   labels = c("a", "b"),
                                   rel_widths = c(1, 2)))

ggsave(plot = FinalFigure, 
       filename = plotfile, 
       width = 7, height = 3)





