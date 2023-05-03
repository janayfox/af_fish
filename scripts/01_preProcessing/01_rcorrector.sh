#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=5-00:00
#SBATCH --account=def-barrett

#######################################################
### Goal: Run rcorrector on raw reads
### Author: Janay Fox
#######################################################

# $1 = R1 file 
# $2 = R2 file 
 
module load rcorrector

perl /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/rcorrector/1.0.4/bin/run_rcorrector.pl -t 12 -1 $1 -2 $2 -od /home/janayfox/scratch/afFishRNA/rcorrector_output
