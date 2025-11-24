#!/usr/bin/env Rscript

### GobySatDNA: exploring the annotation of satellite DNA in the goby genome
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/04/06
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2)
library(cowplot)
library(stringr)
# ============================

# IMPORTANT: set path of the repository
setwd("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Manuscripts/12_Goby/Repositories/GobySD/9_Others")

### Reference fGobFla1
seqlensref <- read.table("data/fGobFla1_ref_seqlen.txt", header = FALSE)
names(seqlensref) <- c("Scaffold", "Size")

# Annotation reference
gff_sdref <- read.table("../../Dryad/Annotation/fGobFla1-GobyTideCluster_v1.00.gff3", header = FALSE)

### TH1 assembly
seqlens <- read.table("data/PomflaTH1xy_seqlen.txt", header = FALSE)
names(seqlens) <- c("Scaffold", "Size")

# Annotation TH1
gff_sd <- read.table("../../Dryad/Annotation/PomflaTH1xy-GobyTideCluster_v1.00.gff3", header = FALSE)

## Output 
barplotsatDNA <- "results/FigS6_BarPlotSatDMA_SDvsRest.png"

# ============================
# Functions
# ============================
## Process the GFF3 into a more palatable table, vectorized version: satDNA
gff2satDNA <- function(gff){
  names(gff) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Extract Name and ID using regex
  gff <- gff %>%
    mutate(
      Name = str_extract(attributes, "Name=[^;]+") %>% stringr::str_replace("Name=", ""),
      ID = str_extract(attributes, "ID=[^;]+") %>% stringr::str_replace("ID=", ""),
      TRC = ifelse(str_detect(Name, "TRC_"),
                   str_match(Name, "(TRC_\\d+)")[,2],
                   NA)) %>% 
    select(Scaffold, start, end, Name, ID, TRC)
  
  return(gff)
}

## Process the GFF3 into a more palatable table, vectorized version: TEs
gff2TEs <- function(gff){
  names(gff) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Extract Name and ID using regex
  gff <- gff %>%
    mutate(
      Name = str_extract(attributes, "Name=[^;]+") %>% str_replace("Name=", ""),
      ID = str_extract(attributes, "ID=[^;]+") %>% str_replace("ID=", "") 
      ) %>% 
    select(Scaffold, start, end, Name, ID)
  return(gff)
}

## Process the GFF3 into a more palatable table, vectorized version: TEs
gff2TEs <- function(gff){
  names(gff) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Extract Name and ID using regex
  gff <- gff %>%
    mutate(
      Name = str_extract(attributes, "Name=[^;]+") %>% str_replace("Name=", ""),
      ID = str_extract(attributes, "ID=[^;]+") %>% str_replace("ID=", "") 
    ) %>% 
    select(Scaffold, start, end, Name, ID)
  return(gff)
}

## Process the GFF3 into a more palatable table, vectorized version: Augustus genes
gff2augustus <- function(gff){
  names(gff) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Extract Name and ID using regex
  gff <- gff %>%
    mutate(
      ID = str_extract(attributes, "ID=[^;]+") %>% str_replace("ID=", ""),
      Parent = str_extract(attributes, "Parent=[^;]+") %>% str_replace("Parent=", "")
    ) %>% 
    select(Scaffold, type, start, end, ID, Parent)
  return(gff)
}

# ============================
# Processing
# ============================
# SD scaffolds
malescf <- "h1tg000060l"
femalescf <- "h2tg000053l" # "ptg000042l" in the primary assembly

# ============================
# Abundance of the most common satDNA in the SD scaffold across the genome
# ============================
gffdfsat <- gff2satDNA(gff_sd)

# Add the sizes of the scaffolds and remove mitochondrion
gffsatDNA <- merge(gffdfsat, seqlens, by = "Scaffold") %>% 
  filter(Scaffold != "PomflaTH1_mt") %>% 
  mutate(range = end - start)

scfcov_sd <- gffsatDNA %>% group_by(Scaffold, Name) %>% summarise(cov = sum(range))
scfcov_sd <- merge(scfcov_sd, seqlens, by = "Scaffold")

# For plotting it's better to call ptg000042l as the hap2 variant
scfcov_sd <- scfcov_sd %>% 
  mutate(Scaffold  = if_else(Scaffold == "ptg000042l", "h2tg000053l", Scaffold))

# Define custom color mapping
my_colors <- c(
  "TRC_1"   = "#FE8776",
  "TRC_3"   = "#3DBEAC",
  "sf10__TRC_75"  = "slateblue3",
  "sf11__TRC_42"  = "#FF00FF",
  "TRC_105" = "#AE123A" )

SDscf_satDNA <- ggplot(
  scfcov_sd %>% filter(Scaffold %in% c(malescf, femalescf), cov > 2000),
  aes(x = Name, y = cov, fill = Name) ) +
  geom_bar(stat = "identity") +
  facet_grid(Scaffold ~ .) +
  scale_fill_manual(
    values = my_colors,
    na.value = "grey80"  # all others not in my_colors will be grey
  ) +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  labs(fill = "Repeat family", title = "Haplotype versions of the SD scaffold (coverage > 2 Kb)") +
  ylab("Coverage (bp)") + xlab("satDNA family")

# In general
rest_satDNA <- ggplot(
  scfcov_sd %>% filter(!Scaffold %in% c(malescf, femalescf)) %>% filter(cov > 50000),
  aes(x = Name, y = cov, fill = Name) ) +
  geom_bar(stat = "identity") +
  # facet_wrap(. ~ Scaffold) +
  scale_fill_manual(
    values = my_colors,
    na.value = "grey80"  # all others not in my_colors will be grey
  ) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  labs(fill = "satDNA family", title = "TH1 genome excluding the SD scaffold (coverage > 50 Kb)") +
  ylab("Coverage (bp)") + xlab("satDNA family")
# So the most common satDNA families in the genome are TRC_1 and TRC_105.


# ============================
# satDNA along the chromosomes
# ============================
gffdfsatref <- gff2satDNA(gff_sdref)

# Add the sizes of the scaffolds and remove mitochondrion
gffsatDNAref <- merge(gffdfsatref, seqlensref, by = "Scaffold") %>% 
  filter(Scaffold != "OZ251433.1") %>% 
  mutate(range = end - start)

ref_to_chr <- c(
  "OZ251410.1" = "1",
  "OZ251411.1" = "2",
  "OZ251412.1" = "3",
  "OZ251415.1" = "6",
  "OZ251413.1" = "4",
  "OZ251414.1" = "5",
  "OZ251416.1" = "7",
  "OZ251417.1" = "8",
  "OZ251418.1" = "9",
  "OZ251419.1" = "10",
  "OZ251420.1" = "11",
  "OZ251421.1" = "12",
  "OZ251422.1" = "13",
  "OZ251423.1" = "14",
  "OZ251424.1" = "15",
  "OZ251425.1" = "16",
  "OZ251426.1" = "17",
  "OZ251427.1" = "18",
  "OZ251428.1" = "19",
  "OZ251429.1" = "20",
  "OZ251430.1" = "21",
  "OZ251431.1" = "22",
  "OZ251432.1" = "23",
  "OZ251433.1" = "mt" )

( satDANchrs <- ggplot() + 
    ## Chromosomes
    geom_rect(data = seqlensref %>% filter(Scaffold != "OZ251433.1"),
              aes(ymin = 1, ymax = Size, xmin = 0.1, xmax = -0.1),
              fill = "gainsboro",
              inherit.aes = FALSE,
              show.legend = FALSE) +
    facet_grid(. ~ Scaffold, labeller = labeller(Scaffold = ref_to_chr), switch = "x") +
    ## repeats
    geom_rect(data = gffdfsatref,
              aes(ymin = start, ymax = end, xmin = 0.1, xmax = -0.1,
                  fill = Name),
              inherit.aes = FALSE,
              alpha = 0.8,
              show.legend = FALSE) +
    scale_fill_manual(
      values = my_colors,
      na.value = "grey40"  # all others not in my_colors will be grey
    ) +
    ## colors
    theme_classic() +
    ylab("Chromosomal position (bp)") +
    theme(axis.text = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 9),
          axis.title.y = element_text(color = "black")) )

satDNAplot <- cowplot::plot_grid(SDscf_satDNA, rest_satDNA, satDANchrs,
                                 nrow = 3,
                                 align = "v",
                                 axis = "tblr",
                                 labels = c("a", "b", "c"),
                                 label_size = 20,
                                 rel_heights = c(1, 0.5, 1))

ggsave(plot = satDNAplot, 
       file = barplotsatDNA, 
       width = 8, height = 12)

