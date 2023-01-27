## Plots trees and colors taxa according to lineage affinities

library(ape)
library(phangorn)
library(phytools)
library(RColorBrewer)
library(tidyverse)
library(ggtree)
library(ggrepel)
library(viridisLite)

#read in data
Tree1 <- read.tree("AllSeqsPops.raxml.support")

# assign the groups
groupInfo <- split(Tree1$tip.label, gsub(".*_", "", Tree1$tip.label))
groupInfo <- lapply(groupInfo, function(x) gsub("_.*", "", x))

# Identify samples by research group
Ceresini <- c("As", "Ds", "Ei534i", "Ecrus326", "Ub0", "Ub3", "Ce535i", "Mr0", "Elcan", "Ce642i")
DelPonte <- c("U1", "U2", "U7", "UFVPY")

# group isolates by host genus
Tree1$tip.label <- gsub("_.*", "", Tree1$tip.label)
Tree1 <- groupOTU(Tree1, groupInfo)

# adjust branch lengths for plotting
Tree1$edge.length <- Tree1$edge.length * 500

# plot the tree
ggtree(Tree1,
       layout = "rectangular", 
       branch.length="branch.length",
       aes(color = group, label = gsub(".*_", "", label))) +
       geom_tiplab(offset = 4, size = 0, linesize = 0.25, align = TRUE, color = "black") +
       scale_color_viridis_d(25, option = "turbo", direction = -1) +
       geom_tiplab(geom = "text", size = 2, linesize = 0, align = TRUE, offset = 4) +
       geom_tippoint(aes(subset=(grepl(paste(Ceresini, collapse = "|"), label, fixed = F) == TRUE), x=x+1), size = 2, position = "identity", color = "blue") +
       geom_tippoint(aes(subset=(grepl(paste(DelPonte, collapse = "|"), label, fixed = F) == TRUE), x=x+0.8), size = 2, shape = "triangle", color = "red", ) +
       guides(col = guide_legend(title = "Host Genus", nrow = 25, width = 5, override.aes = aes(size = 2, label = ""))) +
       geom_treescale(x= 5, y = -3) + geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 80))

# save the tree
ggsave("tree.pdf", dpi = 600, width =10, height =24)

