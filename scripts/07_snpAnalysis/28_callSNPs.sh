#!/bin/bash
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=120G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Construct SuperTranscripts
### Author: Janay Fox
#######################################################

# $1 = R1 file 
# $2 = R2 file 
# $3 = name_output 

module load StdEnv/2020
module load python/3.8.2
module load apptainer/1.1.6
#module load nixpkgs/16.09
module load picard/2.26.3
#module load gatk/3.8
module load samtools/1.17
module load star/2.7.9a

# singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
# trinityrnaseq.v2.15.1.simg /usr/local/bin/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
# --st_fa /lustre04/scratch/janayfox/af_fish_RNA/BA_super_trinity_genes.fasta \
# --st_gtf /lustre04/scratch/janayfox/af_fish_RNA/BA_super_trinity_genes.gtf \
# -p $1 $2 -o $3 \
# --maxram 400000000000 -t 3


singularity exec -e --env-file envfile -B /lustre04/scratch/janayfox/af_fish_RNA \
trinityrnaseq.v2.15.1.simg /lustre04/scratch/janayfox/af_fish_RNA/trinityrnaseq-v2.15.1/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa /lustre04/scratch/janayfox/af_fish_RNA/BN_super_trinity_genes.fasta \
--st_gtf /lustre04/scratch/janayfox/af_fish_RNA/BN_super_trinity_genes.gtf \
--maxram 400000000000 -t 3 \
-samples_file /lustre04/scratch/janayfox/af_fish_RNA/sample_lists/samples_BN_sal_bf.txt \
-o /lustre04/scratch/janayfox/af_fish_RNA/variant_calling/BN

#while read -r value1 value2 value3; do sbatch 27_callSNPs.sh $value1 $value2 $value3; done < BA_files.txt