#!/usr/bin/env Rscript

# Gggenomes_goby: Plotting the sex-determining scaffold of the two-spotted goby *Pomatoschistus flavescens*
#############################################################################
# Sources:
# - https://thackl.github.io/gggenomes/articles/emales.html
# - https://thackl.github.io/gggenomes/reference/gggenomes.html
# - Emil's starfish implementation: https://github.com/egluckthaler/starfish/blob/main/main/locus-viz
# =======================================
# S. Lorena Ament-Velasquez
# 2025-04-07
#############################################################################
# Load the necessary libraries
# ============================
library(tidyverse)
library(gggenomes)
library(ggnewscale)
library(paletteer) # https://r-charts.com/color-palettes/
library(ggthemes)

# ============================
# Read file names
# ============================
# ## Snakemake 
haplotypes <- snakemake@input$fa
gffrepeats <- snakemake@input$satDNA
links <- snakemake@input$paf
gffgenes <- snakemake@input$genes
plotname <- snakemake@output$plot


# # Local
# path2files <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Manuscripts/12_Goby/Repositories/GobySD/1_SDsynteny"
# # Local input
# haplotypes <- paste0(path2files, "/data/haplotypes.fa")
# links <- paste0(path2files, "/minimap/haplotypes_short.paf")
# gffrepeats <- paste0(path2files, "/data/SDscaffolds.gff3")
# # gffgenes <- paste0(path2files, "/data/TH1_PflaAmhr2y.gff3")
# gffgenes <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Manuscripts/12_Goby/Repositories/Dryad/Annotation/TH1_PflaAmhr2y.gff3"
# 
# # # Output files
# plotname <- paste0(path2files, "/results/SDsynteny.pdf")

# ============================
# Parse files
# ============================
# parse sequence length
haplo_seqs <- read_seqs(haplotypes, parse_desc=FALSE) 
haplorepeats <- read_feats(gffrepeats, parser="read_gff3")
suppressWarnings(haplo_links <- read_paf(links, max_tags = 24)) # otherwise Snakemake fails
haplogenes <- read_feats(gffgenes, parser="read_gff3")

# Calculate fraction of identity
haplo_links <- haplo_links %>% mutate(Identity = map_match/map_length)

haplorepeats <- haplorepeats %>%
  mutate(
    family = case_when(
      str_detect(name, "^sf") ~ str_extract(name, "(?<=^sf)[A-Za-z0-9]+(?=__)"),
      str_detect(name, "^TRC_") ~ str_extract(name, "(?<=TRC_)[A-Za-z0-9]+"),
      TRUE ~ NA_character_
    ),
    len = end - start
  )

# Visualize
lowcolor <- "darkslategray1"
highcolor <- "deepskyblue4"

## Change the name of the primary contig to match the haplotype2 (identical) contig
haplo_links <- haplo_links %>%
  mutate(
    seq_id  = if_else(seq_id == "ptg000042l", "h2tg000053l", seq_id),
    seq_id2 = if_else(seq_id2 == "ptg000042l", "h2tg000053l", seq_id2)
  )
haplo_seqs <- haplo_seqs %>% mutate(seq_id  = if_else(seq_id == "ptg000042l", "h2tg000053l", seq_id))
haplorepeats <- haplorepeats %>% mutate(seq_id  = if_else(seq_id == "ptg000042l", "h2tg000053l", seq_id))

## Make a palette for the satDNA families
families <- sort(unique(haplorepeats$family))
family_colors <- paletteer_c("grDevices::Dark 3", length(families))

set.seed(3) # just to make it reproducible
family_colors <- sample(family_colors)
names(family_colors) <- families
# Color specific families
family_colors["1"] <- "#FE8776"
family_colors["3"] <- "#3DBEAC"
family_colors["10"] <- "slateblue3"
family_colors["11"] <- "#FF00FF"
family_colors["105"] <- "#AE123A"

# packageVersion("gggenomes")
# packageVersion("ggplot2")
print(haplo_seqs)

## Plot
syntenyp <- gggenomes(seqs = haplo_seqs,
                      genes = haplogenes %>% filter(type == "CDS"),
                      feats = haplorepeats,
                      links = haplo_links %>% filter(map_length > 5000)) +
  geom_seq(linewidth = 6, alpha = 1, color = "gray80") +  # draw seqs boundaries
  geom_link(aes(colour = Identity, fill = Identity), alpha = 1, size = 0.05, offset = 0.1) + # draw links (offset - how close the links are to the scaffolds)
  geom_gene() +
  scale_fill_gradient(low = lowcolor, high = highcolor) +
  scale_colour_gradient(low = lowcolor, high = highcolor, guide = "none") + # remove the double legend
  ggnewscale::new_scale_fill() + # Geoms added to a plot after this function will use a new scale definition.
  ggnewscale::new_scale_colour() + # Geoms added to a plot after this function will use a new scale definition.
  geom_bin_label() + # add labels of sequences and scale
  geom_feat(aes(colour = family), linewidth = 4, show.legend = FALSE, position = "identity") +  # thickness of TEs
  scale_x_bp(suffix = "b", sep = " ") +
  scale_color_manual(values = family_colors) +
  theme_gggenomes_clean(base_size = 12) +
  theme(legend.title = element_text(size=10), legend.text = element_text(size=8)) + # to make it closer to the base_size
  scale_y_continuous(expand = expansion(mult = 0.1))# Expand the space between scaffolds

ggsave(plot = syntenyp, file = plotname, width = 8, height = 2)

