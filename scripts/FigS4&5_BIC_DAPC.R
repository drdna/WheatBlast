library("adegenet")
library("poppr")
library("RColorBrewer")

n <- 20
palette <- distinctColorPalette(n)

df <- read.table("~/B71v2sh_SNPs_BIC20.gz", row.names = 1, header = FALSE) # StructureIn.str is OK

D <- df2genind(df, ncode = 1, ploidy = 1, NA.char = "-9")

PS <- popsub(D)

PC <- find.clusters(PS, max.n.clust = 50, n.iter = 1e6)

# PC: Chose 80 PCs for this analysis

maxK <- 50

# nrow specifies # iterations for BIC distribution analysis

myMat <- matrix(nrow=50, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(D, n.pca = 200, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

# Plot distribution of K after nrow iterations

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1

# set k range for DAPC

my_k <- 10

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(D, n.pca = 200, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(D, pop = grp_l[[i]]$grp, n.pca = 200, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp

# Plot the cluster memberships (Fig. B)

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data$State
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

#re-order data for tidy plotting

plotOrder <- my_df[my_df$Posterior==1,]
rownames(plotOrder) <- plotOrder$Isolate
plotOrder <- plotOrder[order(plotOrder[,'Isolate']), ]
plotOrder <- plotOrder[order(plotOrder[,'Group']), ]
plotOrder["K"][plotOrder["K"] == 10] <- "A"
strainOrder = plotOrder$Isolate

# gather Ceresini sample data
PCplots <- read.table("~/PlotCeresini", header = F)
colnames(PCplots) <- c("Isolate", "K", "Group", "Posterior")
PCplots["Group"][PCplots["Group"] == 0] <- 11
PCplots["K"][PCplots["K"] == 0] <- "B"
PCplots["K"][PCplots["K"] == 1] <- "B"
PCplots$Group <- factor(PCplots$Group)
rownames(PCplots) <- PCplots$Isolate
PCplots <- PCplots[order(match(rownames(PCplots), strainOrder)), , drop = FALSE]

# gather Farman sample data
MFplots <- read.table("~/PlotFarman", header = F)
colnames(MFplots) <- c("Isolate", "K", "Group", "Posterior")
MFplots["Group"][MFplots["Group"] == 0] <- 11
MFplots["K"][MFplots["K"] == 0] <- "C"
MFplots["K"][MFplots["K"] == 1] <- "C"
MFplots$Group <- factor(MFplots$Group)
rownames(MFplots) <- MFplots$Isolate
MFplots <- MFplots[order(match(rownames(MFplots), strainOrder)), , drop = FALSE]

# merge datasets
plotOrder = rbind(plotOrder, PCplots, MFplots)

# build the plots
facetLabels <- c("All Samples","Near Wheat","Away from Wheat")
names(facetLabels) <- c("A", "B", "C")

p3 <- ggplot(plotOrder, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity", width = 0.5)
p3 <- p3 + scale_x_discrete(guide=guide_axis(angle=90), limits = strainOrder)
p3 <- p3 + facet_grid(K ~ ., space = "fixed", labeller = labeller(K = facetLabels))
p3 <- p3 + ylab("Population membership")
p3 <- p3 + theme(aspect.ratio = 1/20, panel.spacing = unit(1, "pt"), strip.text.y = element_text(size = 12, angle = 90))
p3 <- p3 + scale_fill_manual(values = c("#FF63B6", "#F8766D", "#64B200", "#AEA200", 
                                        "#00A6FF", "#B385FF", "#DB8E00", "#00BD5C",
                                        "#00C1A7", "#00BADE", "black"), breaks = c(1,2,3,4,5,6,7,8,9,10)) 
p3 <- p3 + theme(plot.margin = margin(1,1,1,1, "cm"), panel.border = element_rect(colour = "black", fill = NA),
                 axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                 axis.ticks.x = element_line(size = 0.1), axis.text.x = element_text(angle = 90, hjust = 1, size = 4.5),
                 panel.spacing.y=unit(0.01,"in"), panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(), panel.grid.major.x = element_line(size = 0.1),
                 panel.grid.minor.x = element_blank(), legend.position = "top", )
p3 <- p3 + guides(fill = guide_legend(nrow = 1, byrow = TRUE))
p3

pdf("~/FigS5_DAPC10.pdf", 11, 8.5)
p3
dev.off()

