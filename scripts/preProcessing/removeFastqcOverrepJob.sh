#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett

# $1 = R1 file 
# $2 = R2 file 
# $3 = fastqc file R1 
# $4 = fastqc file R2
 
module load python
python removeFastqcOverrep.py -1 $1 -2 $2 -fql $3 -fqr $4