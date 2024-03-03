#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Run DEG analysis using edgeR in Trinity
### Author: Janay Fox
#######################################################

module load StdEnv/2020
module load apptainer/1.1.6

#run with genes, min CPM, 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon_new/matrix/BN_bf_new_sal.gene.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt --output BN_salmon_DEG_gene_minCPM_new \
# --min_reps_min_cpm '7,1'

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.gene.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt --output BA_salmon_DEG_gene_minCPM \
# --min_reps_min_cpm '7,1' 

#try rerunning with different filter cpm
singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/assemblyReadsBeforeRm/BN/salmon_new/matrix/BN_bf_new_sal.gene.counts.matrix \
--method edgeR --samples_file /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt --output BN_salmon_DEG_gene_minCPM_new \
--min_reps_min_cpm '7,2'

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.gene.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt --output BA_salmon_DEG_gene_minCPM \
# --min_reps_min_cpm '7,2' 