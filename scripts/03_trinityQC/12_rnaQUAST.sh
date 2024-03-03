#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-00:10
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Run RNA QUAST on both assemblies
### Author: Janay Fox
#######################################################

module purge 
module load StdEnv/2020
module load gcc/9.3.0
module load python/3.8.10 
module load quast/5.0.2

python /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/quast/5.0.2/bin/quast.py -o quast.output -t 8 -l "Trinity_trim.fa" -e --rna-finding BN_bf.Trinity.fasta

python /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/quast/5.0.2/bin/quast.py -o quast.output -t 8 -l "Trinity_trim.fa" -e --rna-finding BA_bf.Trinity.fasta