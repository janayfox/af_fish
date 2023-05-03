#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-10:00
#SBATCH --account=def-barrett

#######################################################
### Goal: Trim reads using TrimGalore
### Author: Janay Fox
#######################################################

# $1 = R1 file 
# $2 = R2 file 
 
module purge
module load fastqc
module load python/3.10
source cutAdaptENV/bin/activate

#run for BA
/home/janayfox/scratch/afFishRNA/TrimGalore-0.6.6/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads_BA --length 36 -q 5 --stringency 1 -e 0.1 $1 $2

#run for BN
/home/janayfox/scratch/afFishRNA/TrimGalore-0.6.6/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads_BN --length 36 -q 5 --stringency 1 -e 0.1 $1 $2
