#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Quality check replicates and samples
### Author: Janay Fox
#######################################################

#have to run this way to have access to R packages 
module load singularity/3.8

### going to have to potentially edit sample files here? and make new ones for the uncleaned

#compare replicates 
## reads before removal of overrep ##
#kallisto 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BA/matrix_kallisto_BA_bf/BA_bf_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BN/matrix_kallisto_BN_bf/BN_bf_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# #salmon
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BA/matrix_salmon_BA_bf/BA_bf_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BN/matrix_salmon_BN_bf/BN_bf_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# new BN 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/kallisto/matrix/BN_bf_new_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_kal_bf.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# new BA
singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
--matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/kallisto/matrix/BA_bf_new_kal.isoform.counts.matrix \
--samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal_bf.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt --log2 --CPM --min_rowSums 10 --compare_replicates

#compare across all samples 
#correlation matrix
## reads before removal of overrep ##
#kallisto 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BA/matrix_kallisto_BA_bf/BA_bf_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BN/matrix_kallisto_BN_bf/BN_bf_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_kal.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# #salmon
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BA/matrix_salmon_BA_bf/BA_bf_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BN/matrix_salmon_BN_bf/BN_bf_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

#new BN
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/kallisto/matrix/BN_bf_new_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_kal_bf.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

#new BA
singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
--matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/kallisto/matrix/BA_bf_new_kal.isoform.counts.matrix \
--samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal_bf.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# #PCA
# ## reads before removal of overrep ##
# #kallisto 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BA/matrix_kallisto_BA_bf/BA_bf_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BN/matrix_kallisto_BN_bf/BN_bf_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3

# #salmon
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BA/matrix_salmon_BA_bf/BA_bf_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR  \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BN/matrix_salmon_BN_bf/BN_bf_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3

#new BN
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/kallisto/matrix/BN_bf_new_kal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_kal_bf.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3

#new BA
singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
--matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/kallisto/matrix/BA_bf_new_kal.isoform.counts.matrix \
--samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal_bf.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.counts.matrix \
# --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt --log2 --CPM --min_rowSums 10 --center_rows --prin_comp 3