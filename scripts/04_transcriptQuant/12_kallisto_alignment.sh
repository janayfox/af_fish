#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load singularity/3.8
module load nixpkgs/16.09
module load gcc/9.3.0
module load intel/2020.1.217
module load kallisto/0.46.1

#prepare reference for alignment 
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BA.Trinity.fasta --est_method kallisto \
--trinity_mode --prep_reference

#run alignment and abundance estimation
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BA.Trinity.fasta --seqType fq --SS_lib_type RF --samples_file /data/samples_BA_kal.txt \
--est_method kallisto --trinity_mode --output_dir kallisto_output