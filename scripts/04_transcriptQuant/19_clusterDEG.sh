#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Get preliminary clustering figures
### Author: Janay Fox
#######################################################

module load singularity/3.8

# #BN
# #kallisto
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# -m-atrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/kallisto/matrix/BN_bf_new_kal.isoform.TMM.EXPR.matrix \
# -P 0.001 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_kal_bf.txt

# #salmon
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.TMM.EXPR.matrix \
# -P 0.001 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt

#BA
#kallisto
# singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/kallisto/matrix/BA_bf_new_kal.isoform.TMM.EXPR.matrix \
# -P 0.001 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal_bf.txt

# #salmon
# singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.TMM.EXPR.matrix \
# -P 0.001 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt

#new DEG analysis, focusing on salmon 
# singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.gene.TMM.EXPR.matrix \
# -P 0.001 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt

singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
/lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.gene.TMM.EXPR.matrix \
-P 0.001 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt
