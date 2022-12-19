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
3. Return to terminal #1 and type run to start:
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
1. Reads were aligned to the B71 reference genome and haplotypes were called using the BWT2-GATK.sh script:
```bash
for f in `ls RAW_READS/*fq.gz | awk -F '[/.]' '{print $2}'`; do sbatch $script/BWT2-GATK.sh 70-15.fasta RAW_READS $f; done
```
2. Generate fasta files for MonsterPlex data (reorder based on chromosome, position):
```bash
perl MonsterPlex2Fasta_noMGG.pl 70-15.fasta MPcoverage.bed MPLEX_VCFs > MG_MPLEX.fasta
```
3. Create maximum likelihood tree using RAxML:
```bash
raxml -T 4 -p 48556 -f a -x 48556 -s MonsterPlex_final.fasta -n MonsterPlex_final.raxml -m GTRGAMMA -# 100
```
4. Add support values to nodes


