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

1. Use [Generate_SNP_dataset.pl](/scripts/Generate_SNP_dataset.pl) script to generate .fasta file from SNPs called by iSNPcaller. In this case, the SNPs were called using alignments between each strain and a single reference genome - [B71v2.fasta](https://www.ncbi.nlm.nih.gov/genome/62?genome_assembly_id=1571357):
```bash
perl Generate_SNP_dataset.pl strain_list.txt SNP_COUNTS B71v2sh.fasta B71v2sh_ALIGN_STRINGs > AllSeqs.fasta
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
![FigS2_MLtree.tif](/images/FigS2MLtree.tif)

## Figure S3 & S4. Identifying populations using k-means clustering

1. Generate a down-sampled dataset (every 20th SNP) in STRUCTURE format using [Generate_STRUCTURE.pl](/scripts/Generate_STRUCTURE.pl) script:
```bash
Usage: Generate_STRUCTURE.pl <strain-list> <iSNPcaller-outdir> <ref-genome> <downsample-factor>
```
```bash
Generate_STRUCTURE.pl StrainList B71v2sh_SNPs B71v2sh.fasta 20
awk 'NR > 2' B71v2sh_SNPs_structure20 | gzip > B71v2sh_SNPs_BIC20.gz
```
2. Use B71v2sh_SNPs_BIC20.gz as input to the [FigS3&4_BIC_DAPC.R](/scripts/FigS3&4_BIC_DAPC.R) script:
![FigS3_BIC1-50.png](/images/FigS3_BIC1-50.png)
### FigS3.Population division using discriminant analysis of principal components
![FigS4_DAPC10.tif](/images/FigS4_DAPC10.tif)
### Figure S4. Population memberships of isolates included in this study
 
## Figure S5. Assessing iSNPcaller error rates
Briefly, errors rates were determined by using iSNPcaller to identify SNPs between independent assemblies of the same genome generated from the same raw read dataset.
The [trim-velvet-FR.sh](/scripts/trim-velvet-FR.sh) script was used to perform the following operations:
1. Low quality sequence and adapters were trimmed from the raw reads using Trimmomatic with parameters: ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN 130
2. The resulting paired and unpaired forward reads were concatenated into forwardReads.fq
3. The resulting paired and unpaired reverse reads were concatenated into reverseReads.fq
4. The forward and reverse reads were assembled separately using VelvetOptimiser using the kmer range 89 to 129, with a step size of 2:
```bash
for f in `ls /project/farman_uksr/PE_datasets2/ | awk -F '/|_' '{print $1}' \ 
 | sort | uniq | grep -v GB | grep -v ERR| grep -v ^WB`; \
 do sbatch $script/trim-velvet-FRreads.sh \
 /project/farman_uksr/PE_datasets2 $f yes; done
```
5. Sequence headers were standardized:
```bash
for f in `ls ASSEMBLIES/*fasta`; do perl SimpleFastaHeaders_SB.pl $f; done
```
6. The genomes were repeat-masked using the RMSA_MT.sh SLURM script:
```bash
sbatch RMSA_MT.sh ASSEMBLIES
```
7. Masked genomes generated using the forward and reverse reads were then blasted against one another:
```bash
mkdir SELFBLASTs
cd ASSEMBLIES
for f in `ls *R_nh_masked*`; do blastn -query $f -subject ${f/-R_/-F_} -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../SELFBLASTs/${f/_*/}.${f/R_*/F}.BLAST; done
```
8. SNPs were called using the iSNPcaller SNP calling module:
```bash
perl Run_SU4.pl SELFBLASTs SELFSNPs
```
9. The SNP summary file was condensed into a format suitable for loading into an R script for plotting:
```bash
perl Consense_SNP_summary.pl 
```
10. The condensed SNP summary file was then read into the [FigS5_ClonalvFR_SNPs.R](/scripts/FigS5_ClonalvFR_SNPs.R) script for plotting:

![FigS5_ClonalvFR_SNPs.tif](/images/FigS5_ClonalvFR_SNPs.tif)
### Figure S5. Determining clonality through divergence analysis.




