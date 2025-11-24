#!/usr/bin/env Rscript

# Dotplot of Goby alignments
# =======================================
# Sandra Lorena Ament Velasquez
# 2025-05-28
# =======================================

# ============================
# Load the necessary libraries
# ============================
library(ggplot2, quietly = TRUE)
library(dplyr, warn.conflicts = FALSE)
library(scales) # for the palettes

# ============================
# Input files
# ============================
## Local
# Input
coords <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/IFB/06_Map2Reference/mummer/PomflaTH1hap2-vs-fGobFla1.filter.coords", header = F)
chrlens <- read.table("/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/IFB/06_Map2Reference/data/reference/fGobFla1_ref_seqlen.txt", header = F)
# Output
bigpicturefile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/IFB/06_Map2Reference/results/PomflaTH1hap2-vs-fGobFla1_Dotplot.png"
chr16plotfile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/IFB/06_Map2Reference/results/PomflaTH1hap2-vs-fGobFla1_Dotplot_chr16.png"
chr16barfile <- "/Users/lorena/Library/CloudStorage/Dropbox/VRwork/Analyses/7_Goby/IFB/06_Map2Reference/results/PomflaTH1hap2-vs-fGobFla1_chr16_bar.pdf"


## Snakemake
# Input
coords <-  read.table(snakemake@input$coords, header = F)
chrlens <-  read.table(snakemake@input$chrlens, header = F)
# Output
bigpicturefile <- snakemake@output$big
chr16plotfile <- snakemake@output$chr16
chr16barfile <- snakemake@output$chr16bar

# ============================
# Processing
# ============================
names(coords) <- c("S1", "E1", "S2", "E2", "LEN_R", "LEN_Q", "IDY", "COV_R", "COV_Q", "REF", "QUERY")
names(chrlens) <- c("REF", "LenREF")
  
# Map the chromosome names to the accession numbers in Genbank
ref_to_chr <- c(
  "OZ251410.1" = "chr1",
  "OZ251411.1" = "chr2",
  "OZ251412.1" = "chr3",
  "OZ251415.1" = "chr6",
  "OZ251413.1" = "chr4",
  "OZ251414.1" = "chr5",
  "OZ251416.1" = "chr7",
  "OZ251417.1" = "chr8",
  "OZ251418.1" = "chr9",
  "OZ251419.1" = "chr10",
  "OZ251420.1" = "chr11",
  "OZ251421.1" = "chr12",
  "OZ251422.1" = "chr13",
  "OZ251423.1" = "chr14",
  "OZ251424.1" = "chr15",
  "OZ251425.1" = "chr16",
  "OZ251426.1" = "chr17",
  "OZ251427.1" = "chr18",
  "OZ251428.1" = "chr19",
  "OZ251429.1" = "chr20",
  "OZ251430.1" = "chr21",
  "OZ251431.1" = "chr22",
  "OZ251432.1" = "chr23",
  "OZ251433.1" = "mt" )

coords$Chromosome <- ref_to_chr[coords$REF]
chrlens$Chromosome <- ref_to_chr[chrlens$REF]

# # Assign orientation to the query contigs vs the reference
# coords$orientation <- ifelse(
#   as.numeric(coords$S1 - coords$E1) * as.numeric(coords$S2 - coords$E2) > 0,
#   "forward", "reverse")

coords <- merge(coords, chrlens, by = c("REF", "Chromosome"))

# I'm not interested in the mitochondria
coords <- coords %>% filter(Chromosome != "mt") 

## Sort the chromosome names with natural sorting
chr_levels <- unique(coords$Chromosome)
chr_nums <- as.numeric(gsub("chr", "", chr_levels))

# Get naturally sorted, unique chromosome names
sorted_chr <- chr_levels[order(chr_nums)] # use the chromosome numbers as indexes for ordering
coords$Chromosome <- factor(coords$Chromosome, levels = sorted_chr)

# ============================
# Functions
# ============================

dotplotter2 <- function(data, linesize = 1.5, minaligncov = 100000, minidentity = 96, minLENQ = 10000){
  ## First, select the contigs that have long alignments
  datacov <- data %>% group_by(QUERY, REF) %>% summarise(totalalign = sum(LEN_Q))
  selectedhits <- datacov %>% filter(totalalign > minaligncov) %>% pull(QUERY) # Just take the contigs with decent alignment lengths
  
  data <- merge(data, datacov, by = c("QUERY", "REF")) 
  
  # Asign colors before ordering the contigs so there is more contrast
  n_colors <- length(selectedhits)
  my_colors <- hue_pal()(n_colors)  # Generates `n` evenly spaced colors
  names(my_colors) <- selectedhits       # Name each color by contig
  
  ## Next, order the contigs to match the reference
  # Compute first hit (or earliest reference coordinate) per query contig
  query_order <- data %>%
    filter(QUERY %in% selectedhits, IDY > minidentity, LEN_Q > minLENQ) %>%
    group_by(REF, QUERY) %>%
    summarize(min_ref = min(pmin(S1, E1)), .groups = "drop") %>%
    arrange(REF, min_ref) %>% 
    pull(QUERY) %>% unique() %>% rev()
  
  # Recover the lost contigs
  query_order <- c(setdiff(selectedhits, query_order), query_order )
  
  # Change the factors to match the order
  data$QUERY <- factor(data$QUERY, levels = query_order)
  
  
  dotp <- ggplot(data %>% filter(QUERY %in% selectedhits)) +
    # Alignment segments
    geom_segment(aes(x = S1, xend = E1, y = S2, yend = E2, color = QUERY), linewidth = 1.5) +
    theme_minimal() +
    facet_grid(QUERY ~ REF, scales="free", space = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.8),
          # panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          strip.text.y = element_text(angle = 0), # make strip labels horizontal
          axis.line = element_line(color='black'),
          panel.spacing = unit(0.05, "lines"), # make the facets closer 
          legend.position="none") +
    scale_color_manual(values = my_colors) +
    # scale_color_manual(values = c("forward" = "cyan3", "reverse" = "coral2")) +
    labs(x = "Reference Coordinate (bp)", y = "Query Coordinate (bp)") 
  return(dotp)
}

# ============================
# Scaffolds with sex markers
# ============================
sexscfs <- c("h2tg000004l", "h2tg000009l", "h2tg000023l", "h2tg000034l", "h2tg000038l", "h2tg000053l", "h2tg000074l")
sexscfsals <- coords %>% filter(QUERY %in% sexscfs)

# I want the names of the chromosomes rather than the GenBank accession
sexscfsals$REF2 <- sexscfsals$REF
sexscfsals$REF <- sexscfsals$Chromosome

relevantchrs <- c("chr16", "chr21", "chr3", "chr18")

dotplotter2(sexscfsals %>% filter(Chromosome %in% relevantchrs), 
            minLENQ = 10000, minaligncov = 1000) +
  theme(axis.text.x = element_text(angle = 90) )

  ## Where exactly are the SNPs?
# h2tg000004l are 4 SNPs between 18198203 and 18198362 (159 bp) - chr21
# h2tg000009l are 7 SNPs between 4311506 and 4311673 (167 bp) - chr18
# h2tg000023l are 4 SNPs between 23833541 and 23833694 (153 bp) - chr3
# h2tg000034l are 2 SNPs of coordinates 453815 and 453874 (59 bp) (chr16 anyway) -- contig is 750605 bp
# h2tg000038l are 3 SNPs between 1756287 and 1756442 (155 bp) - chr3
# h2tg000074l is one SNP at 1106662

sexscfsals %>% filter(QUERY == "h2tg000034l", LEN_Q > 5000, S2 > 453815 - 30000, E2 < 453874 + 30000) 
# chr16
sexscfsals %>% filter(QUERY == "h2tg000004l", S2 > 18198203 - 100, E2 < 18198362 + 100) 
# chr21
# Based on that the SNPs should fall somewhere in OZ251430.1:7615330-7615572 -- TRC_47 == Goby1-137, clearly a satellite
sexscfsals %>% filter(QUERY == "h2tg000009l", S2 > 4311506 - 100, E2 < 4311673 + 100) 
# Not a repeat area... chr18
sexscfsals %>% filter(QUERY == "h2tg000023l", S2 > 23833541 - 100, E2 < 23833694 + 100)
# Not a repeat area... chr3
sexscfsals %>% filter(QUERY == "h2tg000038l", S2 > 1756287 - 100, E2 < 1756442 + 100) 
# chr3

# ============================
# Plotting chr16 only
# ============================
### Chr16 in particular (sex chr)
chr16 <- coords %>% filter(Chromosome == "chr16")

# I want the names of the chromosomes rather than the GenBank accession
chr16$REF2 <- chr16$REF
chr16$REF <- chr16$Chromosome

# Remove these because they are distracting
scfsWithSatDNAalignments <- c()

chr16plot <- dotplotter2(chr16 %>% filter(!QUERY %in% scfsWithSatDNAalignments), 
                         minLENQ = 5000) + theme(axis.text.x = element_text(angle = 0))

ggsave(plot = chr16plot + xlab("French goby chr16 (bp)") + ylab("Norwegian goby scaffolds (bp)"), 
       filename = chr16plotfile,
       width = 5, height = 5)

# ============================
### To get a big picture you need to make some sacrifices, so first remove small scfs
# ============================
bigQcoords <- coords %>% filter(COV_Q > 3000000)
chr_nums_bigQ <- as.numeric(gsub("chr", "", bigQcoords$Chromosome))
bigQcoords$REF2 <- bigQcoords$REF
bigQcoords$REF <- chr_nums_bigQ

(bigpicture <- dotplotter2(bigQcoords, minLENQ = 30000, minaligncov = 250000, linesize = 1.1) + 
  theme(axis.text.x = element_blank(),
        # strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        # panel.grid.major = element_blank()
        ) + 
    labs(x = "French goby chromosome", y = "Norwegian goby scaffold") ) 

ggsave(plot = bigpicture, 
       filename = bigpicturefile, 
       width = 7, height = 6)

# ============================
# Where are those translocations exactly?
# ============================

# dotplotter2(bigQcoords %>% filter(Chromosome %in% c("chr3", "chr7", "chr11", "chr20")), 
#             minLENQ = 30000, minaligncov = 250000, linesize = 1.1)

# # ============================
# # Lazy version that is not a proper dotplot
# # ============================
# ggplot(coords) +
#   # Chromosome base line
#   geom_segment(data = distinct(coords, Chromosome, LenREF),
#                aes(x = 0, xend = LenREF, y = -10000, yend = -10000),
#                color = "gray50", size = 1) +
#   
#   # Alignment segments
#   geom_segment(aes(x = S1, xend = E1, y = S2, yend = E2, color = orientation)) +
#   
#   facet_wrap(~ Chromosome) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.8)) +
#   scale_color_manual(values = c("forward" = "cyan3", "reverse" = "coral2")) +
#   labs(x = "Reference Coordinate (bp)", y = "Query Coordinate (bp)")
# 
# 

# 
# ## Create a palette for the different contigd
# contigs <- unique(chr16$QUERY)
# library(scales)
# 
# 
# ggplot(chr16) +
#   # Chromosome base line
#   geom_segment(data = distinct(chr16, Chromosome, LenREF),
#                aes(x = 0, xend = LenREF, y = -10000, yend = -10000),
#                color = "gray50", size = 1) +
#   # Alignment segments
#   geom_segment(aes(x = S1, xend = E1, y = S2, yend = E2, color = QUERY), size = 1.1) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.8), 
#         legend.position="none") +
#   scale_color_manual(values = my_colors) +
#   # scale_color_manual(values = c("forward" = "cyan3", "reverse" = "coral2")) +
#   labs(x = "Reference Coordinate (bp)", y = "Query Coordinate (bp)")
# 
# ## ----- Can I make it more dotplot-like? -------
# 
# chr16cov <- chr16 %>% group_by(QUERY) %>% summarise(totalalign = sum(LEN_Q))
# selectedhits <- chr16cov %>% filter(totalalign > 100000) %>% pull(QUERY)
# 
# #### Order the contigs to match the reference
# # Compute first hit (or earliest reference coordinate) per query contig
# query_order <- chr16 %>%
#   filter(QUERY %in% selectedhits, LEN_R > 10000) %>%
#   group_by(QUERY) %>%
#   summarize(min_ref = min(pmin(S1, E1))) %>%  # handle reverse alignments too
#   arrange(min_ref) %>%
#   pull(QUERY)
# 
# chr16$QUERY <- factor(chr16$QUERY, levels = rev(query_order))
# 
# 
# ggplot(chr16 %>% filter(QUERY %in% selectedhits)) +
#   # Alignment segments
#   geom_segment(aes(x = S1, xend = E1, y = S2, yend = E2, color = QUERY), size = 1.5) +
#   theme_minimal() +
#   facet_grid(QUERY~., scales="free_y", space = "free_y") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.8),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.y = element_blank(),
#         strip.text.y = element_text(angle = 0), # make strip labels horizontal
#         axis.line = element_line(color='black'),
#         legend.position="none") +
#   scale_color_manual(values = my_colors) +
#   # scale_color_manual(values = c("forward" = "cyan3", "reverse" = "coral2")) +
#   labs(x = "Reference Coordinate (bp)", y = "Query Coordinate (bp)") 
# 

