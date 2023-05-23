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
library(ggtext)

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

# set a color palette
c25 <- c(
  "yellow3", # Avena
  "#E31A1C", # red, Bromus
  "green4", # Cynodon
  "#6A3D9A", # purple, Cenchrus
  "#FF7F00", # orange, Digitaria
  "black", # Eleusine
  "gold1", # Echinochloa
  "skyblue2", # Elionurus
  "#FB9A99", # lt pink, Eragrostis
  "palegreen2", # Hakonechloa
  "#FDBF6F", # lt orange, Hordeum
  "#CAB2D6", # lt purple,  Lolium
  "gray70", # Leptochloa
  "khaki2", # Leersia
  "maroon", # Luziola
  "orchid1", # Melinis
  "deeppink1", # Oryza
  "dodgerblue2", # Panicum
  "steelblue4", # Paspalum
  "darkturquoise", # Poa
  "green1", # Setaria
  "yellow4", # Stenotaphrum
  "blue1", # Triticum
  "darkorange4", # Urochloa
  "brown" # Zea
)

# plot the tree
t <- ggtree(Tree1,
       layout = "rectangular", 
       branch.length="branch.length",
       aes(label = gsub(".*_", "", label))) +
  
       # the following block holds key tricks for creating decent legends
       geom_tiplab(aes(color = group), size = 2, linesize = 0.25, align = TRUE, offset = 2, show.legend=F) +
       scale_color_manual(values = c25) +
       geom_polygon(aes(fill = group, x = 0, y = 0)) +
       scale_fill_manual(name = "Host Genus",
                         labels = c("*Avena*",
                                    "*Bromus*",
                                    "*Cenchrus*",
                                    "*Cynodon*",
                                    "*Digitaria*",
                                    "*Eleusine*",
                                    "*Echinochloa*",
                                    "*Elionurus*",
                                    "*Eragrostis*",
                                    "*Hakonechloa*",
                                    "*Hordeum*",
                                    "*Lolium*",
                                    "*Leptochloa*",
                                    "*Luziola*",
                                    "*Leersia*",
                                    "*Melinis*",
                                    "*Oryza*",
                                    "*Panicum*",
                                    "*Paspalum*",
                                    "*Poa*",
                                    "*Setaria*",
                                    "*Stenotaphrum*",
                                    "*Triticum*",
                                    "*Urochloa*",
                                    "*Zea*"),
                         values = c25,
                         guide=guide_legend(nrow=25, keywidth=1,
                                            keyheight=1,
                                            order=1,
                                            override.aes=list(size=5,alpha=1))) +
       # end of block
  
       geom_tippoint(aes(subset=(grepl(paste(Ceresini, collapse = "|"), label, fixed = F) == TRUE), x=x+1), 
                     size = 2, position = "identity", color = "blue") +
       geom_tippoint(aes(subset=(grepl(paste(DelPonte, collapse = "|"), label, fixed = F) == TRUE), x=x+0.8),
                     size = 2, shape = "triangle", color = "red", ) +
       geom_treescale(x= 5, y = -4) + xlim_tree(12) +
       theme(legend.position = c(0.1, 0.88), legend.text = element_markdown())

print(t)

# save the tree
ggsave("FigS3_MLtree.pdf", dpi = 300, width =8, height = 24)
