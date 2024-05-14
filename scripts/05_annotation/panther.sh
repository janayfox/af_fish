#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=100G
#SBATCH --time=5-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Prep for Panther
### Author: Janay Fox
#######################################################

module load StdEnv/2020
module load nixpkgs/16.09
#module load panther/14.1
module load gcc/7.3.0
module load hmmer/3.2.1

./pantherScore2.2.pl -l ../PANTHER18.0_hmmscoring.tgz -D B \
-V -i ../../trinity_output/final_assembly/BA/BA_bf.Trinity.fasta \
-o panther_BA.tsv -n -s

./pantherScore2.2.pl -l ../PANTHER18.0_hmmscoring.tgz -D B \
-V -i ../../trinity_output/final_assembly/BN/BN_bf.Trinity.fasta \
-o panther_BA.tsv -n -s
