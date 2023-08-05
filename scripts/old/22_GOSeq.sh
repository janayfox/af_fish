#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-02:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Run Gene Ontology Analysis 
### Author: Janay Fox
#######################################################

module load StdEnv/2020
module load apptainer/1.1

#extract GO assignments 
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/util/extract_GO_assignments_from_Trinotate_xls.pl \
# --Trinotate_xls /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_trinotate.xls \
# -G --include_ancestral_terms > BN_go_annotations.txt

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/util/extract_GO_assignments_from_Trinotate_xls.pl \
# --Trinotate_xls /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_trinotate.xls \
# -G --include_ancestral_terms > BA_go_annotations.txt

#generate gene lengths file
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/fasta_seq_length.pl \
# /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta > BN_Trinity.fasta.seq_lens

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/TPM_weighted_gene_length.py \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta.gene_trans_map \
# --trans_lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_Trinity.fasta.seq_lens \
# --TPM_matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.isoform.TMM.EXPR.matrix > BN.Trinity.gene_lengths.txt

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/fasta_seq_length.pl \
# /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BA/BA_bf.Trinity.fasta > BA_Trinity.fasta.seq_lens

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/misc/TPM_weighted_gene_length.py \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BA/BA_bf.Trinity.fasta.gene_trans_map \
# --trans_lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_Trinity.fasta.seq_lens \
# --TPM_matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.isoform.TMM.EXPR.matrix > BA.Trinity.gene_lengths.txt

#run GO analysis on DEG 
# singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon/matrix/BN_bf_new_sal.gene.TMM.EXPR.matrix \
# -P 0.01 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BN_sal_bf.txt --examine_GO_enrichment \
# --GO_annots /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
# --gene_lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

# singularity exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/Analysis/DifferentialExpression/analyze_diff_expr.pl \
# --matrix /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/matrix/BA_bf_new_sal.gene.TMM.EXPR.matrix \
# -P 0.01 -C 2 --samples /lustre04/scratch/janayfox/afFishRNA/samples_BA_sal_bf.txt --examine_GO_enrichment \
# --GO_annots /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_go_annotations.txt \
# --gene_lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA.Trinity.gene_lengths.txt

# #run GOseq on gene lists
# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BN_venn_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BN_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BN_k_cluster2_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BN_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BN_30perc_cluster7_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BN_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
/lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
--factor_labeling /lustre04/scratch/janayfox/afFishRNA/BN_mfuzz_cluster5_factor.txt \
--background /lustre04/scratch/janayfox/afFishRNA/BN_background.txt \
--GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
--lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BN_mfuzz_cluster6_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BN_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BN_StF.vs.SwF_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BN_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BN_SwF.vs.SwP_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BN_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate/BN_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate/BN.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BA_venn_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BA_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BA_k_cluster3_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BA_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA.Trinity.gene_lengths.txt

apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
/lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
--factor_labeling /lustre04/scratch/janayfox/afFishRNA/BA_30perc_cluster5_factor.txt \
--background /lustre04/scratch/janayfox/afFishRNA/BA_background.txt \
--GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_go_annotations.txt \
--lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA.Trinity.gene_lengths.txt

apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
/lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
--factor_labeling /lustre04/scratch/janayfox/afFishRNA/BA_mfuzz_cluster5_factor.txt \
--background /lustre04/scratch/janayfox/afFishRNA/BA_background.txt \
--GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_go_annotations.txt \
--lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BA_StF.vs.SwF_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BA_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA.Trinity.gene_lengths.txt

# apptainer exec -e --env-file /lustre04/scratch/janayfox/afFishRNA/envfile -B /lustre04/scratch/janayfox/afFishRNA \
# /lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
# --factor_labeling /lustre04/scratch/janayfox/afFishRNA/BA_SwF.vs.SwP_factor.txt \
# --background /lustre04/scratch/janayfox/afFishRNA/BA_background.txt \
# --GO_assignments /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA_go_annotations.txt \
# --lengths /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA.Trinity.gene_lengths.txt