#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett


module load nixpkgs/16.09
module load StdEnv/2020
module load fastqc/0.11.9

fastqc *.gz