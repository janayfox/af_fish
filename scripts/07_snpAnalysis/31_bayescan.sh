#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=40G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --job-name=BA_filtered
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load StdEnv/2020
module load nixpkgs/16.09 
module load intel/2018.3
module load bayescan/2.1

#filtered 
bayescan_2.1 ./BA_filtered/BA_bayescan_bayescan.txt -od ./BA_filtered -o BA_bayescan_filtered -threads 5 \
-n 10000 -burn 200000 -pr_odds 10000

#bayescan_2.1 ./BN_filtered/BA_bayescan_bayescan.txt -od ./BN_filtered -o BN_bayescan_filtered -threads 5 \
#-n 10000 -burn 200000 -pr_odds 10000

#not filtered for hwe 
#bayescan_2.1 ./BA_full/BA_bayescan_full.txt -od ./BA_full -o BA_bayescan_full -threads 5 \
#-n 10000 -burn 200000 -pr_odds 10000

#bayescan_2.1 ./BN_full/BA_bayescan_full.txt -od ./BN_full -o BN_bayescan_full -threads 5 \
#-n 10000 -burn 200000 -pr_odds 10000