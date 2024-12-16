#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Get Ex90N50 stats on both assemblies
### Author: Janay Fox
#######################################################

## can only run this after transcript quantification - use isoform expression matrix##

module load singularity/3.8

## Get Ex90N50 and Ex90 ##
singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
/data/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.TMM.EXPR.matrix \
/data/readsBeforeRmoverrep/BN/BN_bf.Trinity.fasta transcript | tee BN_bf_new_salmon.ExN50.transcripts.stats

singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
/data/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.TMM.EXPR.matrix \
/data/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta transcript | tee BA_bf_new_salmon.ExN50.transcripts.stats

# ## plot results ## 
singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
/data/BN_bf_new_salmon.ExN50.transcripts.stats

singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
/data/BA_bf_new_salmon.ExN50.transcripts.stats

# estimate TPM thresholds for transcript counting and filtering ## 
singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
--E_inputs /data/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.TMM.EXPR.matrix.E-inputs 

singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
--E_inputs /data/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.TMM.EXPR.matrix.E-inputs 
