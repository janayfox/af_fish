#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

## can only run this after transcript quantification ##

module load singularity/3.8

singularity exec -e -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/contig_ExN50_statistic.pl \
/data/quant_output_cleanedReads/kallisto_output/BA/matrix_kallisto_BA/kallisto_BA.isoform.TMM.EXPR.matrix \
/data/cleanedReads/BA/BA.Trinity.fasta transcript | tee BA_kallisto.ExN50.transcripts.stats

singularity exec -e -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/contig_ExN50_statistic.pl \
/data/quant_output_cleanedReads/kallisto_output/BN/matrix_kallisto_BN/kallisto_BN.isoform.TMM.EXPR.matrix \
/data/cleanedReads/BN/BN.Trinity.fasta transcript | tee BN_kallisto.ExN50.transcripts.stats