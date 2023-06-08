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

module load singularity/3.8

# #kallisto
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BA/matrix_kallisto_BA_bf/BA_bf_kal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal.txt --output BA_kallisto_DEG

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BN/matrix_kallisto_BN_bf/BN_bf_kal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BN_kal.txt --output BN_kallisto_DEG

# #salmon
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BA/matrix_salmon_BA_bf/BA_bf_sal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal.txt --output BA_salmon_DEG

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BN/matrix_salmon_BN_bf/BN_bf_sal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal.txt --output BN_salmon_DEG

#new BN
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/kallisto/matrix/BN_bf_new_kal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BN_kal_bf.txt --output BN_kallisto_DEG

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt --output BN_salmon_DEG

#new BA
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/kallisto/matrix/BA_bf_new_kal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BA_kal_bf.txt --output BA_kallisto_DEG

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt --output BA_salmon_DEG

# re run with genes, min CPM, focussing on salmon
singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon_new/matrix/BN_bf_new_sal.gene.counts.matrix \
--method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt --output BN_salmon_DEG_gene_minCPM_new \
--min_reps_min_cpm '7,1'

singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.gene.counts.matrix \
--method edgeR --samples_file /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt --output BA_salmon_DEG_gene_minCPM \
--min_reps_min_cpm '7,1' 