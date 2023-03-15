#!/bin/bash

#SBATCH --time 48:00:00
#SBATCH --job-name=trim-velvet
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --partition=normal
#SBATCH --mem=500GB
#SBATCH --mail-type ALL
#SBATCH -A coa_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST

dir=$1

f=$2

mkdir $f 

cp $dir/$f*_1*f*q* $f/

cp $dir/$f*_2*f*q* $f/

cd $f

if [ $3 == 'yes' ]
then
  singularity run --app trimmomatic039 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf trimmomatic PE \
  -threads 16 -phred33 -trimlog ${f}_errorlog.txt \
  $f*_1*.f*q* $f*_2*.f*q* \
  ${f}_R1_paired.fq ${f}_R1_unpaired.fq \
  ${f}_R2_paired.fq ${f}_R2_unpaired.fq \
  ILLUMINACLIP:/project/farman_uksr/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:20:20 MINLEN:130;
fi


# hard code read filenames

cat ${f}_R1_paired.fq ${f}_R1_unpaired.fq > forwardReads.fq

cat ${f}_R2_paired.fq ${f}_R2_unpaired.fq > reverseReads.fq

# remove previous failed run data

rm -r velvet_${f}_R

rm -r velvet_${f}_F

# now run velvet on reverse reads

singularity run --app perlvelvetoptimiser226 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf VelvetOptimiser.pl \
 -s 89 -e 129 -x 2 -d velvet_${f}_R -f ' -short -fastq reverseReads.fq'

mv velvet_${f}_R/contigs.fa velvet_${f}_R/${f}"-R.fasta"
cp velvet_${f}_R/${f}"-R.fasta" /scratch/farman/ASSEMBLIES/${f}"-R.fasta"

prefix=`ls velvet_${f}_R/*Logfile.txt`

mv $prefix velvet_${f}_R/${f}_R_${prefix/*\//}

# now run velvet on forward reads

singularity run --app perlvelvetoptimiser226 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf VelvetOptimiser.pl \
 -s 89 -e 129 -x 2 -d velvet_${f}_F -f ' -short -fastq forwardReads.fq'

mv velvet_${f}_F/contigs.fa velvet_${f}_F/${f}"-F.fasta"
cp velvet_${f}_F/${f}"-F.fasta" /scratch/farman/ASSEMBLIES/${f}"-F.fasta"

prefix=`ls velvet_${f}_F/*Logfile.txt`

mv $prefix /scratch/farman/ASSEMBLIES/${f}_F_${prefix/*\//}

# remove temp unneeded files/folders
rm forwardReads.fq
rm reverseReads.fq
rm -r velvet_${f}_F
rm -r velvet_${f}_R
rm ../*gz
rm ../*fq
rm -r ../auto*
cd ..
