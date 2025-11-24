#!/usr/bin/env Rscript

### PixyGobySD: Divergence along the SD scaffold between males and females
#############################################################################

# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/04/26
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2)
library(stringr) # stringr::str_match
library(ggpubr)
library(cowplot)
# ============================

# ============================
# Input files
# ============================
## Snakemake
# Annotation
gffTE <- read.table(snakemake@input$TEs)
gffsatDNA <- read.table(snakemake@input$satDNA)

# Pi
hetzygosity <- read.table(snakemake@input$het, header = TRUE)
# Dxy
dxysexes <- read.table(snakemake@input$dxy, header = TRUE)

# Output
compoundplot <- snakemake@output$het
dxyplot <- snakemake@output$dxy

# ============================
# Input parameters
# ============================

## Palette
femalecol <- "#FFD247"
malecol <- "#A482C4"

# start of the central satDNA block
midsatDNA <- 31104200
lensexchr <- 34401611

# Sex chromosome
sexchr <- "OZ251425.1"

# Min ammount of genotyped sites
mingenotyped <- 10000
# ============================
# Functions
# ============================
## Process the GFF3 into a more palatable table, vectorized version: satDNA
gff2genenames <- function(gff){
  names(gff) <- c("chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Extract Name and ID using regex
  gff <- gff %>%
    mutate(
      Name = str_extract(attributes, "Name=[^;]+") %>% str_replace("Name=", ""),
      ID = str_extract(attributes, "ID=[^;]+") %>% str_replace("ID=", ""),
      TRC = ifelse(str_detect(Name, "TRC_"),
                   str_match(Name, "(TRC_\\d+)")[,2],
                   NA)) %>% 
    select(chromosome, start, end, Name, ID, TRC)
  
  return(gff)
}

# Function to merge overlapping intervals
merge_intervals <- function(df) {
  merged <- list()
  current_start <- df$start[1]
  current_end <- df$end[1]
  
  for (i in 2:nrow(df)) {
    this_start <- df$start[i]
    this_end <- df$end[i]
    
    if (this_start <= current_end) {
      # Overlaps -> extend the current region
      current_end <- max(current_end, this_end)
    } else {
      # No overlap -> save current and start new
      merged[[length(merged) + 1]] <- c(current_start, current_end)
      current_start <- this_start
      current_end <- this_end
    }
  }
  # Add the last interval
  merged[[length(merged) + 1]] <- c(current_start, current_end)
  
  # Return as a data frame
  do.call(rbind, merged) %>%
    as.data.frame() %>%
    rename(start = V1, end = V2)
}

# ============================
# Processing
# ============================
# Map the chromosome names to the accession numbers in Genbank
ref_to_chr <- c(
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
  "OZ251433.1" = "mt" )

## satDNA annotation
gffsatDNAdf <- gff2genenames(gffsatDNA) # %>% filter(chromosome == "h1tg000060l")


## Prepare the heterozygosity df
hetzygoclean <- hetzygosity %>% mutate(Midpos = ((window_pos_2 - window_pos_1)/2) + window_pos_1) %>%
  mutate(Sex = case_when(pop == "PflaKGBDf" ~ "Female",
                         pop == "PflaHELAf" ~ "Female",
                         pop == "PflaKGBHm" ~ "Male",
                         pop == "PflaHELEm" ~ "Male")) %>% 
  mutate(het = count_diffs/no_sites)

# Remove windows with too few sites
hetzygoclean$het[which(hetzygoclean$no_sites <= mingenotyped )] <- NA

## Prepare the dxy df
dxysexesclean <- dxysexes %>% mutate(Midpos = ((window_pos_2 - window_pos_1)/2) + window_pos_1) %>%
  mutate(het = count_diffs/no_sites)

# Remove windows with too few sites
dxysexesclean$het[which(dxysexesclean$no_sites <= mingenotyped )] <- NA

# ============================
# Pi
# ============================
## ---- Prepare the representatives of the rest of the genome ----
## satDNA
# When printing the pdf, the satDNA get exagerated, so just take the main ones
relevantsatDNA <- c("TRC_1", "TRC_3", "sf10_", "sf11_", "TRC_105")

focalscfs <- hetzygoclean$chromosome %>% unique() %>% sort()
gffsatDNAdf_scfs <- gffsatDNAdf %>% filter(chromosome %in% focalscfs) %>% filter(TRC %in% relevantsatDNA)

# Average heterozygosity (general)
meanhet <- hetzygoclean %>% group_by(pop) %>% summarise(mean_het = mean(het, na.rm = T), sdcov = sd(het, na.rm = T) )

# no real difference
hetzygoclean %>% filter(chromosome != sexchr) %>% group_by(pop) %>% summarise(mean_het = mean(het, na.rm = T), sdcov = sd(het, na.rm = T) )

### ---- Plotting ----
# To plot the SD region
sdregion <- data.frame(chromosome = sexchr, start = midsatDNA, end = lensexchr)
# Proportion of the non-SD region
100 - ((sdregion$end - sdregion$start)*100)/sdregion$end # 90.41495

# Data frame with label info
sdlabel <- data.frame(
  chromosome = "OZ251425.1",
  x = (midsatDNA + lensexchr) / 2,  # midpoint of SD region
  y = 0.018,                      # y position above the black bar
  label = "X region")

# Define custom color mapping
my_colors <- c(
  "TRC_1"   = "#FE8776",
  "TRC_3"   = "#3DBEAC",
  "sf10__TRC_75"  = "#AE123A",
  "sf11__TRC_42"  = "#FF00FF",
  "TRC_105" = "#F8AF50" )


hetscan <- ggplot(hetzygoclean, 
                   aes(x = Midpos, y = het, colour = Sex, shape = pop)) +
    geom_hline(yintercept= mean(meanhet$mean_het), colour = "gray70") +
    geom_point(alpha = 0.5, size = 1) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    theme_bw() +
    facet_grid(chromosome ~ ., labeller = labeller(chromosome = ref_to_chr)) +
    ggtitle("Heterozygosity along representative autosomes and the sex chromosome (Chr16)") +
    theme(strip.background = element_rect(fill = "white"),
          plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Heterozygosity") +
    xlab("Window position") +
    scale_shape(name = "Individual") +
    scale_colour_manual(values= c("Female" = femalecol, "Male" = malecol)) +
    coord_cartesian(ylim = c(0, 0.02)) +
    # Mark the SD region
    geom_rect(data = sdregion,
              aes(xmin = start, xmax = end, ymin = 0.0157, ymax = 0.016),
              colour = "black", fill = "black",
              inherit.aes = FALSE,
              show.legend = FALSE) +
    # Add the label only for OZ251425.1
    geom_text(data = sdlabel,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = 3,
              colour = "black") +
    # Add the satDNA annotation
    geom_rect(data = gffsatDNAdf_scfs,
              aes(xmin = start, xmax = end, ymin = 0, ymax = -0.001, fill = Name),
              inherit.aes = FALSE,
              alpha = 0.3, # 0.5 for png
              show.legend = FALSE) +
    scale_fill_manual(
      values = my_colors,
      na.value = "gray70"  # all others not in my_colors will be grey
    )

## --- Distribution of heterozygosity between the sexes ----

# Prepare significance values based on the 2 pairwise comparisons
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf)/2, symbols = c("****", "***", "**", "*", "ns"))

my_comparisons <- list( c("PflaKGBDf", "PflaKGBHm"), 
                        # c("PflaKGBDf", "PflaHELEm"), 
                        c("PflaKGBDf", "PflaHELAf"), 
                        c("PflaKGBHm", "PflaHELEm"), 
                        c("PflaHELEm", "PflaHELAf"))
hetzygoclean %>% dplyr::summarize(kruskal_p = kruskal.test(avg_pi ~ pop)$p.value)

# ## All chromosomes
# ggplot(hetzygoclean, # %>% filter(Midpos < midsatDNA), 
#        aes(y = het, fill = Sex, x = pop)) +
#   geom_hline(yintercept= mean(meanhet$mean_het), colour = "gray70") +
#   geom_violin() +
#   geom_jitter(width = 0.1, alpha = 0.3) +
#   facet_grid(chromosome ~ .) +
#   scale_y_continuous(limits = c(0, 0.025)) +
#   theme_classic() + ylab("Heterozygosity") + xlab("Individual") +
#   theme(legend.position="none", axis.title.x = element_blank()) +
#   scale_fill_manual(values= c("Female" = femalecol, "Male" = malecol)) +
#   stat_compare_means(comparisons = my_comparisons, 
#                      p.adjust.method = "bonferroni", 
#                      method = "wilcox.test",
#                      label = "p.adj", 
#                      label.y = c(0.016, 0.018, 0.021, 0.023),
#                      symnum.args = symnum.args)


# Contrasting the chr16 SD region
(hetSDviolin <- ggplot(hetzygoclean %>% filter(Midpos > midsatDNA, chromosome == sexchr), 
                       aes(y = het, fill = Sex, x = pop)) +
    geom_hline(yintercept= mean(meanhet$mean_het), colour = "gray70") +
    geom_violin() +
    geom_jitter(width = 0.1, alpha = 0.3) +
    scale_y_continuous(limits = c(0, 0.025)) +
    ggtitle("Heterozygosity in the X region") +
    theme_classic() + ylab("Heterozygosity") + xlab("Individual") +
    theme(legend.position="none", axis.title.x = element_blank(),
          plot.title = element_text(size=11, hjust = 0.5, face = "bold")) +
    scale_fill_manual(values= c("Female" = femalecol, "Male" = malecol)) +
    stat_compare_means(comparisons = my_comparisons, 
                       p.adjust.method = "bonferroni", 
                       method = "wilcox.test",
                       label = "p.adj", 
                       label.y = c(0.016, 0.018, 0.021, 0.023),
                       symnum.args = symnum.args)) # Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent

(hetNoSDviolin <- ggplot(hetzygoclean %>% filter(Midpos < midsatDNA, chromosome == sexchr), 
                         aes(y = het, fill = Sex, x = pop)) +
    geom_hline(yintercept= mean(meanhet$mean_het), colour = "gray70") +
    geom_violin() +
    geom_jitter(width = 0.1, alpha = 0.3) +
    scale_y_continuous(limits = c(0, 0.025)) +
    ggtitle("Heterozygosity in the rest of Chr16") +
    theme_classic() + ylab("Average heterozygosity") + 
    ylab(" ") + # I want the empty space so that the plot label doesn't overalp with the axis
    theme(legend.position="none", axis.title.x = element_blank(),
          plot.title = element_text(size=11, hjust = 0.5, face = "bold")) +
    scale_fill_manual(values= c("Female" = femalecol, "Male" = malecol)) +
    stat_compare_means(comparisons = my_comparisons,
                       p.adjust.method = "bonferroni", 
                       method = "wilcox.test",
                       label = "p.adj", 
                       label.y = c(0.016, 0.018, 0.021, 0.023),
                       symnum.args = symnum.args) )# Add significance levels # I can't tell if it's doing the same tests, but the plotting is consistent


hetviolins <- cowplot::plot_grid(hetSDviolin, hetNoSDviolin, 
                                 labels = c('b', 'c'), 
                                 label_size = 15, 
                                 align = "tblr")

(compoundp <- cowplot::plot_grid(hetscan, hetviolins, 
                        labels = c('a', NA), 
                        rel_heights = c(1,0.6),
                        label_size = 15, 
                        align = "tblr", nrow = 2))

ggsave(plot = compoundp, filename = compoundplot, width = 9, height = 7)

# ============================
# Compare heterozygosity between males and females (It's not super clear so I didn't use this plot in the end)
# ============================

het_diff <- hetzygoclean %>% #filter(pop %in% c("PflaKGBHm", "PflaKGBDf")) %>% 
  group_by(chromosome, Midpos) %>%
  summarise(
    male_het = mean(het[Sex == "Male"], na.rm = TRUE),
    female_het = mean(het[Sex == "Female"], na.rm = TRUE)
  ) %>%
  mutate(diff = male_het - female_het)

ggplot(het_diff, 
       aes(x = Midpos, y = diff)) +
  geom_hline(yintercept= 0, colour = "gray40") +
  geom_point(alpha = 0.5, size = 1) +
  geom_line(alpha = 0.8, linewidth = 0.8) +
  theme_bw() +
  facet_grid(chromosome ~ .) +
  ylab("Heterozygosity difference\nbetween males and females") +
  xlab("Window position") +
  scale_shape(name = "Individual") +
  
  # Mark the SD region
  geom_rect(data = sdregion,
            aes(xmin = start, xmax = end, ymin = -0.0012, ymax = -0.0014),
            colour = "red", fill = "red",
            inherit.aes = FALSE,
            show.legend = FALSE) +
  # Add the satDNA annotation
  geom_rect(data = gffsatDNAdf_scfs,
            aes(xmin = start, xmax = end, ymin = -0.01, ymax = -0.011, fill = TRC),
            inherit.aes = FALSE,
            alpha = 0.5,
            show.legend = FALSE)


# ============================
# Dxy
# ============================

(dxyscan <- ggplot(dxysexesclean, 
                   aes(x = Midpos, y = het)) +
   geom_point(alpha = 0.5, size = 1) +
   geom_line(alpha = 0.8, linewidth = 0.8) +
   theme_bw() +
   facet_grid(chromosome ~ .) +
   ylab("Heterozygosity") +
   xlab("Window position") +
   scale_shape(name = "Individual") +
   scale_colour_manual(values= c("Female" = femalecol, "Male" = malecol)) +
   # Mark the SD region
   geom_rect(data = sdregion,
             aes(xmin = start, xmax = end, ymin = -0.0012, ymax = -0.0014),
             colour = "black", fill = "black",
             inherit.aes = FALSE,
             show.legend = FALSE) +
   # annotate("text", x=midsatDNA/2, y= -0.0017, label= "SD region", colour = "black") +
   # Add the satDNA annotation
   geom_rect(data = gffsatDNAdf_scfs,
             aes(xmin = start, xmax = end, ymin = 0, ymax = -0.001, fill = TRC),
             inherit.aes = FALSE,
             alpha = 0.5,
             show.legend = FALSE))

# There is no differencein dxy... it's like just mirroring the heterozygosity
ggsave(plot = dxyscan, filename = dxyplot, width = 5, height = 5)

