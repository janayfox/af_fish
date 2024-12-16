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

module load StdEnv/2020
module load apptainer/1.1.6

# singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon_new/matrix/BN_bf_new_sal.gene.TMM.EXPR.matrix \
# -P 0.01 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt

# singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.gene.TMM.EXPR.matrix \
# -P 0.01 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt

#try new 
# singularity exec -e --env-file /lustre04/scratch/janayfox/af_fish_RNA/envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.gene.TMM.EXPR.matrix \
# -P 0.01 -C 2 --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt

singularity exec -e --env-file /lustre04/scratch/janayfox/af_fish_RNA/envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
/lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.gene.TMM.EXPR.matrix \
-P 0.01 -C 2 --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt