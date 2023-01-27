# WheatBlast
Data and scripts for the manuscript: A Re-evaluation of Phylogenomic Data Reveals that Current Understanding in Wheat Blast Population Biology and Epidemiology is Obfuscated by Oversights in Population Sampling

## Figure 1. Generation of a distance tree
1. SNPs were called in all x all pairwise fashion using iSNPcaller.
```bash
perl iSNPcaller.pl Poryzae
```
2. Open new terminal window and copy genomes into iSNPcaller directory:
```bash
cp RAW_GENOMES/*fasta Poryzae/GENOMES
```
3. Return to terminal #1 and type run to start. Briefly, iSNPcaller will perform the following operations:
  i)   Convert sequence header to standard format.
  ii)  Repeat mask all genomes
  iii) Run BLAST to align all genomes in pairwise fashion
  iv)  Call SNPs only in uniquely aligned regions of the genome
```bash
run
```
4. Generate pairwise distance matric in MEGA format:
```bash
perl Pairwise_matrix_MEGA.pl Poryzae/SNP_COUNTS/SNP_counts_XXXXXXXXXX > SNPcounts_WheatBlast.meg
```
5. Use MEGA X to generate a neighbor-joining tree, rendered in radiation format. Edit in Adobe Illustrator.

![Poryzae_distance_tree.tiff](/data/Poryzae_distance_tree.tiff)

## Figure S1. Generation of sampling map

  

## Figure S2. Generation of a maximum likelihood tree

1. Use  .pl script to generate .fasta file from SNPs called by iSNPcaller.pl

