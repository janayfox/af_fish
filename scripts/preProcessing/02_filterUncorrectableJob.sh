#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-02:00
#SBATCH --account=def-barrett

# $1 = R1 file 
# $2 = R2 file 
# $3 = sample ID 
 
module load python
python filterUncorrectable.py -1 $1 -2 $2 -s $3
