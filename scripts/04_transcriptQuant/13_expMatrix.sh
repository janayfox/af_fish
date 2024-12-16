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

## Create matrices ## 
singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
/lustre04/scratch/janayfox/afFishRNA/trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
--est_method salmon --name_sample_by_basedir --out_prefix BN_bf_new_sal \
--gene_trans_map /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta.gene_trans_map \
--quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BN/salmon_new/BN_bf_salquant_files.txt

singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/afFishRNA \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
--est_method salmon --name_sample_by_basedir --out_prefix BA_bf_new_sal \
--gene_trans_map /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta.gene_trans_map \
--quant_files /lustre04/scratch/janayfox/afFishRNA/quant_output/assemblyReadsBeforeRm/BA/salmon/BA_bf_salquant_files.txt