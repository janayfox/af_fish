#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-01:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Generate NxExN50 Stats for both assemblies
### Author: Janay Fox
#######################################################

module load singularity/3.8

#for BA
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/TrinityStats.pl /data/BA.Trinity.fasta

singularity exec -e -B /home/janayfox/scratch/afFishRNA/readsBeforeRmoverrep/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/TrinityStats.pl /data/BA_bf.Trinity.fasta

# #for BN
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BN:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/TrinityStats.pl /data/BN.Trinity.fasta

singularity exec -e -B /home/janayfox/scratch/afFishRNA/readsBeforeRmoverrep/BN:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/TrinityStats.pl /data/BN_bf.Trinity.fasta