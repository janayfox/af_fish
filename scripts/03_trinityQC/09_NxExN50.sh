#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-01:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load singularity/3.8

#for BA
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/TrinityStats.pl /data/BA.Trinity.fasta

#for BN
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BN:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/TrinityStats.pl /data/BN.Trinity.fasta