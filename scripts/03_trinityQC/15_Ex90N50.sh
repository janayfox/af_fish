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

# ## Get Ex90N50 and Ex90 ##
# cleaned reads
# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_cleanedReads/kallisto_output/BA/matrix_kallisto_BA/BA_cl_kal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BA/BA.Trinity.fasta transcript | tee BA_cl_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_cleanedReads/kallisto_output/BN/matrix_kallisto_BN/BN_cl_kal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BN/BN.Trinity.fasta transcript | tee BN_cl_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_cleanedReads/salmon_output/BA/matrix_salmon_BA/BA_cl_sal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BA/BA.Trinity.fasta transcript | tee BA_cl_salmon.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_cleanedReads/salmon_output/BN/matrix_salmon_BN/BN_cl_sal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BN/BN.Trinity.fasta transcript | tee BN_cl_salmon.ExN50.transcripts.stats

# ## pre-cleaning reads mapped to assembly made from cleaned reads 
# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_bfRemOvrRep/kallisto_output/BA/matrix_kallisto_BA_bf/BA_bf_kal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BA/BA.Trinity.fasta transcript | tee BA_bf_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_bfRemOvrRep/kallisto_output/BN/matrix_kallisto_BN_bf/BN_bf_kal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BN/BN.Trinity.fasta transcript | tee BN_bf_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_bfRemOvrRep/salmon_output/BA/matrix_salmon_BA_bf/BA_bf_sal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BA/BA.Trinity.fasta transcript | tee BA_bf_salmon.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/quant_output_bfRemOvrRep/salmon_output/BN/matrix_salmon_BN_bf/BN_bf_sal.isoform.TMM.EXPR.matrix \
# /data/cleanedReads/BN/BN.Trinity.fasta transcript | tee BN_bf_salmon.ExN50.transcripts.stats

# # new BN alignment 
# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/assemblyReadsBeforeRm/BN/kallisto/matrix/BN_bf_new_kal.isoform.TMM.EXPR.matrix \
# /data/readsBeforeRmoverrep/BN/BN_bf.Trinity.fasta transcript | tee BN_bf_new_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
# /data/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.TMM.EXPR.matrix \
# /data/readsBeforeRmoverrep/BN/BN_bf.Trinity.fasta transcript | tee BN_bf_new_salmon.ExN50.transcripts.stats

# new BA alignment 
singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
/data/quant_output/assemblyReadsBeforeRm/BA/kallisto/matrix/BA_bf_new_kal.isoform.TMM.EXPR.matrix \
/data/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta transcript | tee BA_bf_new_kallisto.ExN50.transcripts.stats

singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
/data/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.TMM.EXPR.matrix \
/data/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta transcript | tee BA_bf_new_salmon.ExN50.transcripts.stats


# ## plot results ## 
# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BA_cl_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BN_cl_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BA_cl_salmon.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BN_cl_salmon.ExN50.transcripts.stats


# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BA_bf_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BN_bf_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BA_bf_salmon.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BN_bf_salmon.ExN50.transcripts.stats

# #new BN
# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BN_bf_new_kallisto.ExN50.transcripts.stats

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
# /data/BN_bf_new_salmon.ExN50.transcripts.stats

#new BA
singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
/data/BA_bf_new_kallisto.ExN50.transcripts.stats

singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
/data/BA_bf_new_salmon.ExN50.transcripts.stats

## estimate TPM thresholds for transcript counting and filtering ## 
# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_cleanedReads/kallisto_output/BA/matrix_kallisto_BA/BA_cl_kal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_cleanedReads/kallisto_output/BN/matrix_kallisto_BN/BN_cl_kal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_cleanedReads/salmon_output/BA/matrix_salmon_BA/BA_cl_sal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_cleanedReads/salmon_output/BN/matrix_salmon_BN/BN_cl_sal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_bfRemOvrRep/kallisto_output/BA/matrix_kallisto_BA_bf/BA_bf_kal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_bfRemOvrRep/kallisto_output/BN/matrix_kallisto_BN_bf/BN_bf_kal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_bfRemOvrRep/salmon_output/BA/matrix_salmon_BA_bf/BA_bf_sal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/quant_output_bfRemOvrRep/salmon_output/BN/matrix_salmon_BN_bf/BN_bf_sal.isoform.TMM.EXPR.matrix.E-inputs 

# #new BN
# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.TMM.EXPR.matrix.E-inputs 

# singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
# --E_inputs /data/quant_output/assemblyReadsBeforeRm/BN/kallisto/matrix/BN_bf_new_kal.isoform.TMM.EXPR.matrix.E-inputs 

#new BA
singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
--E_inputs /data/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.TMM.EXPR.matrix.E-inputs 

singularity exec -e --env-file envfile -B /home/janayfox/scratch/afFishRNA:/data \
trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/try_estimate_TPM_filtering_threshold.Rscript \
--E_inputs /data/quant_output/assemblyReadsBeforeRm/BA/kallisto/matrix/BA_bf_new_kal.isoform.TMM.EXPR.matrix.E-inputs 