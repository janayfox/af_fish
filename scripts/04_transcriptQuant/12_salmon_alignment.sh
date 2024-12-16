#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=3-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Use Salmon to quantify transcrtipts
### Author: Janay Fox
#######################################################

module load singularity/3.8
module load gcc/9.3.0
module load StdEnv/2020
moduel load openmpi/4.0.3
module load salmon/1.4.0

# #prepare reference for alignment 
singularity exec -e -B /home/janayfox/scratch/afFishRNA/readsBeforeRmoverrep/BN:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BN_bf.Trinity.fasta --est_method salmon \
--trinity_mode --prep_reference

singularity exec -e -B /home/janayfox/scratch/afFishRNA/readsBeforeRmoverrep/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BA_bf.Trinity.fasta --est_method salmon \
--trinity_mode --prep_reference

# #run alignment and abundance estimation 
singularity exec -e -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/readsBeforeRmoverrep/BN/BN_bf.Trinity.fasta --seqType fq --SS_lib_type RF \
--samples_file /data/readsBeforeRmoverrep/BN/samples_BN_sal_bf.txt \
--est_method salmon --trinity_mode --output_dir salmon_output

singularity exec -e -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta --seqType fq --SS_lib_type RF \
--samples_file /data/readsBeforeRmoverrep/BA/samples_BA_sal_bf.txt \
--est_method salmon --trinity_mode --output_dir salmon_output

#rerun
singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta --seqType fq --SS_lib_type RF \
--samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt \
--est_method salmon --trinity_mode --output_dir salmon_output