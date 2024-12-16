#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-05:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Get preliminary clustering figures
### Author: Janay Fox
#######################################################

module load StdEnv/2020
module load apptainer/1.1.6

# apptainer exec -e -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/util/filter_low_expr_transcripts.pl \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/assemblyReadsBeforeRm/BN/salmon_new/matrix/BN_bf_new_sal.isoform.TMM.EXPR.matrix \
# --transcripts /lustre04/scratch/janayfox/af_fish_RNA/trinity_output/final_assembly/BN/BN_bf.Trinity.fasta \
# --min_expr_any 2.2 | tee BN_bf.filtered.Trinity.fasta

# apptainer exec -e -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/filter_low_expr_transcripts.pl \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.TMM.EXPR.matrix \
# --transcripts /lustre04/scratch/janayfox/af_fish_RNA/trinity_output/final_assembly/BA/BA_bf.Trinity.fasta \
# --min_expr_any 1.1 | tee BA_bf.filtered.Trinity.fasta

# #then I need to generate a new gene-trans-mapping file
# singularity exec -e --env-file /lustre04/scratch/janayfox/af_fish_RNA/envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
# BN_bf.filtered.Trinity.fasta > filtered.BN_bf.Trinity.fasta.gene-trans-map.txt

# singularity exec -e --env-file /lustre04/scratch/janayfox/af_fish_RNA/envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
# BA_bf.filtered.Trinity.fasta > filtered.BA_bf.Trinity.fasta.gene-trans-map.txt

#then redo transcript and gene quant to get new count matrixes
# # #prepare reference for alignment 
# singularity exec -e -B /lustre04/scratch/janayfox/af_fish_RNA\
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
# --transcripts /lustre04/scratch/janayfox/af_fish_RNA/BN_bf.filtered.Trinity.fasta --est_method salmon \
# --trinity_mode --prep_reference

# singularity exec -e -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
# --transcripts /lustre04/scratch/janayfox/af_fish_RNA/BA_bf.filtered.Trinity.fasta --est_method salmon \
# --trinity_mode --prep_reference

# # #run alignment and abundance estimation 
# singularity exec -e -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
# --transcripts /lustre04/scratch/janayfox/af_fish_RNA/BN_bf.filtered.Trinity.fasta --seqType fq --SS_lib_type RF \
# --samples_file /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt \
# --est_method salmon --trinity_mode --output_dir BN_salmon_output_filtered

# singularity exec -e -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
# --transcripts /lustre04/scratch/janayfox/af_fish_RNA/BA_bf.filtered.Trinity.fasta --seqType fq --SS_lib_type RF \
# --samples_file /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt \
# --est_method salmon --trinity_mode --output_dir BA_salmon_output_filtered

# # find . -maxdepth 2 -name "quant.sf" | tee BN_filtered_salquant_files.txt

# ## Create matrices ## 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method salmon --name_sample_by_basedir --out_prefix BN_bf_filtered_sal \
# --gene_trans_map /lustre04/scratch/janayfox/af_fish_RNA/BN_bf.filtered.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BN_filtered/BN_filtered_salquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method salmon --name_sample_by_basedir --out_prefix BA_bf_filtered_sal \
# --gene_trans_map /lustre04/scratch/janayfox/af_fish_RNA/BA_bf.filtered.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BA_filtered/BA_filtered_salquant_files.txt

## Quality check ##
#compare replicates 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BN_filtered/matrix/BN_bf_filtered_sal.gene.counts.matrix \
# --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt --log2 --CPM --min_rowSums 10 --compare_replicates

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BA_filtered/matrix/BA_bf_filtered_sal.gene.counts.matrix \
# --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt --log2 --CPM --min_rowSums 10 --compare_replicates

#compare across all samples 
#correlation matrix
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BN_filtered/matrix/BN_bf_filtered_sal.gene.counts.matrix \
# --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BA_filtered/matrix/BA_bf_filtered_sal.gene.counts.matrix \
# --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt --log2 --CPM --min_rowSums 10 --sample_cor_matrix

# #PCA
# ## reads before removal of overrep ##
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BN_filtered/matrix/BN_bf_filtered_sal.gene.counts.matrix \
# --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt --log2 --CPM --center_rows --prin_comp 4

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/PtR \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BA_filtered/matrix/BA_bf_filtered_sal.gene.counts.matrix \
# --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt --log2 --CPM --center_rows --prin_comp 4

#run DEG analysis with genes, min CPM, 
singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BN_filtered/matrix/BN_bf_filtered_sal.gene.counts.matrix --method edgeR \
--samples_file /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt --output BN_salmon_DEG_gene_minCPM_filtered \
--min_reps_min_cpm '7,1'

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BA_filtered/matrix/BA_bf_filtered_sal.gene.counts.matrix \
# --method edgeR --samples_file /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt --output BA_salmon_DEG_gene_minCPM_filtered \
# --min_reps_min_cpm '7,1' 

# ## cluster ##
# singularity exec -e --env-file /lustre04/scratch/janayfox/af_fish_RNA/envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BN_filtered/matrix/BN_bf_filtered_sal.gene.TMM.EXPR.matrix \
# -P 0.01 -C 2 --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt

# singularity exec -e --env-file /lustre04/scratch/janayfox/af_fish_RNA/envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/af_fish_RNA/quant_output/BA_filtered/matrix/BA_bf_filtered_sal.gene.TMM.EXPR.matrix \
# -P 0.01 -C 2 --samples /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BA_sal_bf.txt