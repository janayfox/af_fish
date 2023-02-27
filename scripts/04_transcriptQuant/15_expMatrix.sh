#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#have to run this way to have access to R packages 
module load StdEnv/2020
module load trinity/2.14.0
module load r/4.2.1

#create matrix
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl \
--est_method kallisto --gene_trans_map ./cleanedReads/BA/BA.Trinity.fasta.gene_trans_map \
--name_sample_by_basedir --quant_files ./BA_kallisto_abundance_files.txt

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl \
--est_method kallisto --gene_trans_map ./cleanedReads/BN/BN.Trinity.fasta.gene_trans_map \
--name_sample_by_basedir --quant_files ./BN_kallisto_abundance_files.txt