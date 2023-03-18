## Plots trees and colors taxa according to lineage affinities

library(ape)
library(phangorn)
library(phytools)
library(RColorBrewer)
library(tidyverse)
library(ggtree)
library(ggrepel)
library(viridisLite)
library(data.table)

#read in data
Tree1 <- read.tree("~/AllSeqs4.fasta.raxml.support")

# assign the groups
groupInfo <- split(Tree1$tip.label, gsub(".*_", "", Tree1$tip.label))

# code to write out IDs of strains used
#labels <- Tree1$tip.label
#write.csv(list(labels), "~/MLTreeStrains.csv")

groupInfo <- lapply(groupInfo, function(x) gsub("_.*", "", x))

# Identify samples by research group
Ceresini <- c("As", "Ds", "Ei534i", "Ecrus326", "Ub0", "Ub3", "Ce535i", "Mr0", "Elcan", "Ce642i")
DelPonte <- c("U1", "U2", "U7", "UFVPY")

# group isolates by host genus
Tree1$tip.label <- gsub("_.*", "", Tree1$tip.label)
Tree1 <- groupOTU(Tree1, groupInfo)

# adjust branch lengths for plotting
Tree1$edge.length <- Tree1$edge.length * 500

#set a color palette
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# plot the tree
ggtree(Tree1,
       layout = "rectangular", 
       branch.length="branch.length",
       aes(label = gsub(".*_", "", label))) +
       geom_tiplab(aes(color = group), size = 2, linesize = 0.25, align = TRUE, offset = 2) +
       scale_color_manual(values = c25) +
       geom_tippoint(aes(subset=(grepl(paste(Ceresini, collapse = "|"), label, fixed = F) == TRUE), x=x+1), 
                     size = 2, position = "identity", color = "blue") +
       geom_tippoint(aes(subset=(grepl(paste(DelPonte, collapse = "|"), label, fixed = F) == TRUE), x=x+0.8),
                     size = 2, shape = "triangle", color = "red", ) +
       guides(col = guide_legend(title = "Host Genus", nrow = 25, width = 2, 
                     override.aes = aes(size = 6, label = ""))) +
       geom_treescale(x= 5, y = -3) + xlim_tree(12) +
       theme(legend.position = c(0.2, 0.85))

# save the tree
ggsave("FigS2_MLtree.pdf", dpi = 300, width =8, height = 24)

