#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett

# $1 = full path to rRNA database (do not include fasta suffix)
# $2 = R1 file 
# $3 = R2 file 
# $4 = sample ID

module load nixpkgs/16.09
module load gcc/7.3.0
module load intel/2018.3
module load bowtie2/2.3.4.3

bowtie2 --nofw --quiet --very-sensitive-local --phred33  -x $1 -1 $2 -2 $3 --threads 12 --met-file ${4}_bowtie2_metrics.txt --al-conc-gz blacklist_paired_aligned_${4}.fq.gz --un-conc-gz blacklist_paired_unaligned_${4}.fq.gz  --al-gz blacklist_unpaired_aligned_${4}.fq.gz --un-gz blacklist_unpaired_unaligned_${4}.fq.gz
