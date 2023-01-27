## Figure 1. Generation of a distance tree
1. SNPs were called in all-by-all pairwise fashion using [iSNPcaller](https://github.com/drdna/iSNPcaller).
```bash
perl iSNPcaller.pl Poryzae
```
2. Open new terminal window and copy genomes into iSNPcaller directory:
```bash
cp RAW_GENOMES/*fasta Poryzae/GENOMES
```
3. Return to terminal #1 and type run to start. Briefly, iSNPcaller will perform the following operations:
  i)   Convert sequence header to standard format;
  ii)  Repeat mask all genomes;
  iii) Run BLAST to align all genomes in pairwise fashion; and
  iv)  Call SNPs only in uniquely aligned regions of the genome
```bash
run
```
4. Generate pairwise distance matric in MEGA format:
```bash
perl Pairwise_matrix_MEGA.pl Poryzae/SNP_COUNTS/SNP_counts_XXXXXXXXXX > SNPcounts_WheatBlast.meg
```
5. Use MEGA X to generate a neighbor-joining tree, rendered in radiation format. Edited in Adobe Illustrator to create figure seen in publication.

![Poryzae_distance_tree.tiff](/data/Poryzae_distance_tree.tiff)

## Figure S1. Generation of sampling map

1. Use [FigS1_SamplingMap.Rmd](/scripts/FigS1_SamplingMap.Rmd) script to build map from [SamplingLocations.xlsx](/data/SamplingLocations.xlsx) file.

![FigS1_SamplingMap.png](/data/FigS1_SamplingMap.png)

## Figure S2. Generation of a maximum likelihood tree

1. Use [Create_fasta_dataset2022.pl](/scripts/Create_fasta_dataset2022.pl) script to generate .fasta file from SNPs called by iSNPcaller. In this case, the SNPs were called using alignments between each strain and a single reference genome - [B71v2.fasta](https://www.ncbi.nlm.nih.gov/genome/62?genome_assembly_id=1571357):
```bash
perl Create_fasta_dataset.pl SNP_COUNTS B71v2sh.fasta B71v2sh_ALIGN_STRINGs > AllSeqs.fasta
```
2. Add host information to each sequence header using [AddPopInfo2Fasta.pl](/scripts/AddPopInfo2Fasta.pl):
```bash
perl AddPopInfo2Fasta.pl StrainHostList.txt AllSeqs.fasta > AllSeqsPops.fasta
```
4. Run RAxML on AllSeqsPops.fasta file (note: fasta file size = 2.2 Gb):
```bash
raxmlHPC-AVX -p 1234 -f a -x 1234 -s AllSeqsPops.fasta -n AllSeqsPops.raxml -m GTRGAMMA #- 1000
```
5. Add support values to nodes:
```bash
raxml -T 2 -f b -m GTRGAMMA -n support -t RAxML_bestTree.AllSeqsPops.raxml -z RAxML_bootstrap.AllSeqsPops.raxml
```
7. Use the resulting .support file to build tree using [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) as implemented in the [FigS2_MLtree.R](/scripts/FigS2_MLtree.R) script:




