#!/usr/bin/env Rscript

### MaleKmers: Looking for SD candidates
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/03/28
# ++++++++++++++++++++++++++++++++++++++++++++++


# ============================
# Load the necessary libraries
# ============================
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2)
# ============================

# ============================
# Input files
# ============================
## Snakemake 
kmercoverage <- read.table(snakemake@input$kmers, header = FALSE)
minN <- snakemake@params$minN

# Output
pointsfile <- snakemake@output$points
chrpaintfile <- snakemake@output$paint

# ============================
# Processing
# ============================
names(kmercoverage) <- c("Scaffold", "Start", "End", "Coverage")

covcounts <- kmercoverage %>% group_by(Scaffold) %>% dplyr::summarise(n = n(), meancov = mean(Coverage) ) #, density = meancov/n

kmercoverage <- kmercoverage %>% mutate(chrnum = factor(Scaffold)) # For chrs painting

# ============================
# Plot male coverage per window along the bigger scaffolds
# ============================

bigbois <- covcounts %>% filter(n > minN) %>% .$Scaffold

## Plot the coverage of male k-mers per window as a scatter plot
(bigbiopoints <- ggplot(kmercoverage %>% filter(Scaffold %in% bigbois), aes(x = End, y = Coverage, colour = Coverage)) +
  geom_hline(yintercept=0.05, linetype='dashed', colour = "gray") +
  geom_point(alpha = 0.6) + 
  facet_wrap(.~Scaffold) +
  theme_classic() +
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Scaffold position"))

ggsave(file = pointsfile, plot = bigbiopoints, width = 8, height = 7)

# ============================
# Chromosome painting version of the plot
# ============================

(chrpaint <- ggplot(data = kmercoverage %>% 
         filter(Scaffold %in% bigbois) %>% 
         mutate(chrnum = factor(chrnum, levels = rev(unique(chrnum)))) ) +
  geom_segment(aes(x = chrnum, 
                   xend = chrnum, 
                   y = Start, yend = End, color = Coverage), linewidth = 3) +
  coord_flip() +
  ylab("Chromosome position") + xlab("Scaffold") +
  scale_color_gradient(low = "khaki1", high = "red") + 
  theme_bw() )

ggsave(file = chrpaintfile, plot = chrpaint, width = 6, height = 10)
