#!/usr/bin/env Rscript

### SexLocusVariants: Plotting the sex-associated variants to our candidate sex-locus 
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/03/28 - 2025/11/04
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2)
library(paletteer) # https://r-charts.com/color-palettes/
library(ggthemes)
# ============================
# Set working directory
setwd("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Manuscripts/12_Goby/Repositories/GobySD/4_PlotMaleKmers")

## HAP1
kmercoverage1 <- read.table("../3_MaleKmers/reports/male_coverage_hap1.txt", header = FALSE)
radsex1 <- read.table("data/radsex_hap1_inv.tsv", header = T)
sexingmarkers <- read.table("data/hap1_SexMarkers_selected_corrected.gff3", header = FALSE )

## HAP2
kmercoverage2 <- read.table("../3_MaleKmers/reports/male_coverage_hap2.txt", header = FALSE)
radsex2 <- read.table("data/radsex_hap2_inv.tsv", header = T)

# Annotation of satDNA on the SD scaffolds
gff <- read.table("../2_SDsynteny/data/SDscaffolds.gff3", header = FALSE)
gffamhr2y <- read.table("data/TH1_PflaAmhr2y.gff3") %>% filter(V3 == "exon")

## Output
hap1vshap2kmerplot <- "results/Fig4_hap1vshap2kmerplot.png"

# ============================
# Functions
# ============================

preparesnps_stacks <- function(df, pval = FALSE){
  if (pval) {
    cleandf <- df %>% dplyr::select(CHR, P)
  } else {
    cleandf <- df %>% dplyr::select(CHR)
  }
  cleandf$POS <- df$BP
  # cleandf$POS <- as.numeric(sub(":.*", "", cleandf$SNP))
  return(cleandf)
}


## Process the GFF3 into a more palatable table
gff2genenames <- function(gff){
  # version 3
  names(gff) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
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
  gff <-  gff %>% select(Scaffold, start, end) %>% cbind(., "Name" = genenames, "ID" = IDs, "TRC" = TRClist) # , "type" = feattypes
  return(gff)
}

# ============================
# Processing
# ============================
# Lengths of the scaffolds
chromosome_lengths <- data.frame(Scaffold = c("h1tg000060l", "h2tg000053l"),
                                 length = c(3958625, 5037559))

# Put them together
kmercoverage <- rbind(kmercoverage1, kmercoverage2) %>% filter(V1 %in% chromosome_lengths$Scaffold)
names(kmercoverage) <- c("Scaffold", "Start", "End", "Coverage")

# Prepare annotation of amhr2y
names(gffamhr2y) <- c("Scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# An upper maximum kmer coverage for plotting
MAXkmerCov <- max(kmercoverage$Coverage)

# Put the radsex markers together
radsex <- rbind(radsex1, radsex2 %>% filter(Contig == "h2tg000053l"))
names(radsex)[names(radsex) == 'Contig'] <- 'Scaffold'
# ============================
# Prepare the annotation
# ============================
## satDNA annotation
gffdf <- gff2genenames(gff) 

# Markers used for sexing
sexymarkers <- gff2genenames(sexingmarkers) 

## Change the name of the primary contig to match the haplotype2 (identical) contig
gffdf <- gffdf %>% mutate(Scaffold = if_else(Scaffold == "ptg000042l", "h2tg000053l", Scaffold), 
                          len = end - start)
sexymarkers <- sexymarkers %>% mutate(Scaffold = if_else(Scaffold == "ptg000042l", "h2tg000053l", Scaffold))

## Make a palette for the satDNA families
families <- sort(unique(gffdf$TRC))
family_colors <- paletteer_c("grDevices::Dark 3", length(families))

set.seed(3) # just to make it reproducible
family_colors <- sample(family_colors)
names(family_colors) <- families
# Color specific families
family_colors["TRC_1"] <- "#FE8776"
family_colors["TRC_3"] <- "#3DBEAC"
family_colors["sf10_"] <- "slateblue3"
family_colors["sf11_"] <- "#FF00FF"
family_colors["TRC_105"] <- "#AE123A"

# But it overwhelms the pdf, snif :(. So only use the main families for plotting
relevantsatDNA <- c("TRC_1", "TRC_3", "sf10_", "sf11_", "TRC_105")

# ============================
# Male-specific kmers coverage + satDNA + sex-associated SNPs
# ============================

# Make a small dataframe for the satDNA label
satDNAlabels <- data.frame(Scaffold = chromosome_lengths$Scaffold,
                           x = c(4200000, 5300000), y = -0.005, label = "satDNA")

# Small dataframe to mark the SD region
SDregion <- data.frame(Scaffold = chromosome_lengths$Scaffold,
                       start = c(1521581, 1698399), 
                       end = chromosome_lengths$length, 
                       label = c("SD region", "Gametologue region")) %>% 
  mutate(midpoint = start + (end - start)/2)

( hap1hap2 <- ggplot(kmercoverage, aes(x = End, y = Coverage, colour = Coverage)) +
    
    ## Scaffold
    geom_rect(data = chromosome_lengths,
              aes(xmin = 1, xmax = length, ymin = -0.004, ymax = -0.001),
              fill = "gray60",
              inherit.aes = FALSE, show.legend = FALSE) +
    facet_wrap(Scaffold ~ ., strip.position = "top", nrow = 2,
               labeller = as_labeller(c("h1tg000060l" = "Haplotype 1 (h1tg000060l)", "h2tg000053l" = "Haplotype 2 (h2tg000053l)"))) +
    geom_point(alpha = 0.9, size = 2) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
          # legend.position = "bottom",
          legend.key.size = unit(0.4, 'cm'),
          legend.title = element_text(size=9),
          legend.position = c(0.9, 0.8),  # (x, y) coordinates inside plotting area
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.title.x=element_blank()) +
    scale_colour_continuous(high = "darkslategrey", low = "darkslategray3") +
    ylab("Male k-mer coverage") +
    ylim(-0.031, MAXkmerCov + 0.025) +

    ## Add the SD region mark
    geom_rect(data = SDregion,
              aes(xmin = start, xmax = end, ymin = 0.15, ymax = 0.151),
              fill = "black", colour = "black",
              inherit.aes = FALSE,
              alpha = 0.5,
              show.legend = FALSE)  +
    # label of SD region mark
    geom_text(data = SDregion,
              aes(x = midpoint, y = 0.16, label = label),
              colour = "black", size = 3, inherit.aes = FALSE) +

    ## Add the satDNA annotation
    geom_rect(data = gffdf %>% filter(TRC %in% relevantsatDNA),
              aes(xmin = start, xmax = end, ymin = -0.010, ymax = -0.005,
                  fill = TRC),
              inherit.aes = FALSE,
              alpha = 0.5, # 0.5 if a png, but lower for pdf
              show.legend = FALSE)  +
    scale_fill_manual(values = family_colors) +
    # label of satDNA
    geom_text(data = satDNAlabels,
              aes(x = x, y = y, label = label),
              colour = "hotpink", size = 3, inherit.aes = FALSE) +
    
    ## amhr2y
    geom_rect(data = gffamhr2y,
              aes(xmin = start, xmax = end, ymin = -0.026, ymax = -0.031),
              colour = "black",
              inherit.aes = FALSE,
              alpha = 0.5,
              show.legend = FALSE)  +
    # Label of amhr2y
    geom_text(
      data = data.frame(
        Scaffold = unique(gffamhr2y$Scaffold),
        x = 4180000, y = -0.028, label = "amhr2y"),
      aes(x = x, y = y, label = label),
      colour = "black", size = 3, fontface = 'italic', inherit.aes = FALSE ) +
    
    ## Add markers used for sexing eggs and larvae
    geom_rect(data = sexymarkers,
              aes(xmin = start, xmax = end, ymin = -0.026, ymax = -0.031),
              colour = "red",
              inherit.aes = FALSE,
              alpha = 0.5,
              show.legend = FALSE)  +
    # Label for sexing markers
    geom_text(data = data.frame(Scaffold = unique(sexymarkers$Scaffold),
        x = 700000, y = -0.028, label = "Markers used for sexing"),
      aes(x = x, y = y, label = label),  
      colour = "red", size = 3, inherit.aes = FALSE ) +
    
    ## Add the SNPs significantly associated to sex (RADSex)
    geom_point(data = radsex %>% mutate(Coverage = 0), 
               aes(x = as.numeric(Position), 
                   y = -0.018), 
               size = 2, alpha = 0.5, colour = "darkslategrey", shape = 2) +
    annotate("text", x=700000, y=-0.016, label= "RADSex markers", colour = "darkslategrey", size = 3) 
)

ggsave(plot = hap1hap2, 
       filename = hap1vshap2kmerplot, 
       width = 6, height = 5)

