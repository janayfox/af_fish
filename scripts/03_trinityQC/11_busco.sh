#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-10:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load StdEnv/2020
module load gcc/9.3.0
module load openmpi/4.0.3
module load busco/5.2.2
module load hmmer/3.3.2

# busco -i ./BA.Trinity.fasta -l ./vertebrata_odb10 -o cleaned_busco_BA -c 4 -m transcriptome --offline

# busco -i ./BN.Trinity.fasta -l ./vertebrata_odb10 -o cleaned_busco_BA -c 4 -m transcriptome --offline

busco -i ./BN_bf.Trinity.fasta -l ../../vertebrata_odb10 -o cleaned_busco_BN_bf -c 4 -m transcriptome --offline