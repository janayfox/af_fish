#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=40G
#SBATCH --time=2-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Construct SuperTranscripts
### Author: Janay Fox
#######################################################

module load StdEnv/2020
module load python/3.8.2
module load scipy-stack/2020a
module load apptainer/1.1.6

apptainer exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
--trinity_fasta ./trinity_output/final_assembly/BA/BA_bf.Trinity.fasta

singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
--trinity_fasta ./trinity_output/final_assembly/BN/BN_bf.Trinity.fasta