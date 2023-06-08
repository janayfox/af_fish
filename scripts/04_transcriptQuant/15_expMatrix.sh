#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-08:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Construct expression matrices
### Author: Janay Fox
#######################################################

module load singularity/3.8

# find . -maxdepth 2 -name "quant.sf" | tee BA_bf_salquant_files.txt
# find . -maxdepth 2 -name "abundance.tsv" | tee BA_bf_kalquant_files.txt

## Create matrices ## 
# # For cleaned read results 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method kallisto --name_sample_by_basedir --out_prefix BA_cl_kal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BA/BA.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_cleanedReads/kallisto_output/BA/BA_cl_kalquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method kallisto --name_sample_by_basedir --out_prefix BA_cl_kal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BA/BA.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_cleanedReads/kallisto_output/BA/BA_cl_kalquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method kallisto --name_sample_by_basedir --out_prefix BN_cl_kal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BN/BN.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_cleanedReads/kallisto_output/BN/BN_cl_kalquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method salmon --name_sample_by_basedir --out_prefix BA_cl_sal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BA/BA.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_cleanedReads/salmon_output/BA/BA_cl_salquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method salmon --name_sample_by_basedir --out_prefix BN_cl_sal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BN/BN.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_cleanedReads/salmon_output/BN/BN_cl_salquant_files.txt

# # For pre-cleaning reads mapped to assembly based on cleaned reads
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method kallisto --name_sample_by_basedir --out_prefix BA_bf_kal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BA/BA.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BA/BA_bf_kalquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method kallisto --name_sample_by_basedir --out_prefix BN_bf_kal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BN/BN.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/kallisto_output/BN/BN_bf_kalquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method salmon --name_sample_by_basedir --out_prefix BA_bf_sal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BA/BA.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BA/BA_bf_salquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method salmon --name_sample_by_basedir --out_prefix BN_bf_sal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/cleanedReads/BN/BN.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/quant_output_bfRemOvrRep/salmon_output/BN/BN_bf_salquant_files.txt

# # For assembly based on reads before cleaning 
# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method kallisto --name_sample_by_basedir --out_prefix BN_bf_new_kal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BN/BN_bf.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/kallisto/BN_bf_kalquant_files.txt

singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
/lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
--est_method salmon --name_sample_by_basedir --out_prefix BN_bf_new_sal \
--gene_trans_map /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta.gene_trans_map \
--quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon_new/BN_bf_salquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method kallisto --name_sample_by_basedir --out_prefix BA_bf_new_kal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/kallisto/BA_bf_kalquant_files.txt

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
# trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
# --est_method salmon --name_sample_by_basedir --out_prefix BA_bf_new_sal \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta.gene_trans_map \
# --quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/BA_bf_salquant_files.txt