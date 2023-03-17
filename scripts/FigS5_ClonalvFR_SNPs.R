#

library(ggplot2)
library(reshape2)
library("PopGenome")

my_df <- read.table("~/Clonal_v_FR_genome_SNPs.txt", header=FALSE)
colnames(my_df) <- c("blast", "SNPs", "type")

pdf("~/FigS5_ClonalvFR_SNPs.pdf", 3, 6)
ggplot(my_df, aes(x= type, y = SNPs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(color = my_df$type, size = 3, width = 0.2, alpha = 0.5) +
  theme_bw() + 
  scale_y_continuous(breaks = seq(0, 90, by = 10)) +
  ylab("SNPs/Mb uniquely\naligned sequence") +
  scale_x_discrete(labels = c("\"clonal\"\nisolates", "false SNP\nerror rate")) +
  theme(axis.text.y =element_text(size=10),
        axis.title.y=element_text(size=12, vjust = 2),
        axis.title.x = element_blank(), axis.text.x=element_text(size=12))
dev.off()

