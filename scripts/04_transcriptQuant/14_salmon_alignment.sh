#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load singularity/3.8
module load gcc/9.3.0
module load StdEnv/2020
moduel load openmpi/4.0.3
module load salmon/1.4.0

#prepare reference for alignment 
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BA.Trinity.fasta --est_method salmon \
--trinity_mode --prep_reference

singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BN:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BN.Trinity.fasta --est_method salmon \
--trinity_mode --prep_reference

#run alignment and abundance estimation on cleaned reads
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BA.Trinity.fasta --seqType fq --SS_lib_type RF --samples_file /data/samples_BA_sal.txt \
--est_method salmon --trinity_mode --output_dir salmon_output

singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BN:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts /data/BN.Trinity.fasta --seqType fq --SS_lib_type RF --samples_file /data/samples_BN_sal.txt \
--est_method salmon --trinity_mode --output_dir salmon_output

#run alignment and abundance estimation on uncleaned reads
singularity exec -e -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/align_and_estimate_abundance.pl \
--transcripts /data/cleanedReads/BA/BA.Trinity.fasta --seqType fq --SS_lib_type RF \
--samples_file /data/readsBeforeRmoverrep/BA/samples_BA_sal_bf.txt \
--est_method salmon --trinity_mode --output_dir salmon_output

singularity exec -e -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/align_and_estimate_abundance.pl \
--transcripts /data/cleanedReads/BN/BN.Trinity.fasta --seqType fq --SS_lib_type RF \
--samples_file /data/readsBeforeRmoverrep/BN/samples_BN_sal_bf.txt \
--est_method salmon --trinity_mode --output_dir salmon_output