#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-05:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Run Trinotate on assemblies
### Author: Janay Fox
#######################################################

module load StdEnv/2020
module load nixpkgs/16.09
module load tmhmm/2.0c
module load python/3.11.2 
module load singularity/3.8

# #for running interactively
# singularity shell -e -B /home/janayfox/scratch:/data trinotate.v4.0.0.simg 
# export TRINOTATE_HOME=/usr/local/src/Trinotate

# #make database 
# singularity exec -e -B /home/janayfox/scratch:/data trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --create --db myTrinotate.sqlite --trinotate_data_dir ./trinotateDatDir --use_diamond

#initialize sqlite database beluga
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate/myTrinotate.sqlite --init \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta.gene_trans_map \
# --transcript_fasta /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta \
# --transdecoder_pep /lustre04/scratch/janayfox/afFishRNA/trinotate/BN/BN_bf.Trinity.fasta.transdecoder.pep

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/myTrinotate.sqlite --init \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BA/BA_bf.Trinity.fasta.gene_trans_map \
# --transcript_fasta /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BA/BA_bf.Trinity.fasta \
# --transdecoder_pep /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA/BA_bf.Trinity.fasta.transdecoder.pep

# #run sequence analyses
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate/myTrinotate.sqlite \
# --transcript_fasta /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta \
# --transdecoder_pep /lustre04/scratch/janayfox/afFishRNA/trinotate/BN/BN_bf.Trinity.fasta.transdecoder.pep \
# --trinotate_data_dir /lustre04/scratch/janayfox/afFishRNA/trinotate/trinotateDatDir \
# --run "swissprot_blastp swissprot_blastx pfam infernal EggnogMapper" --use_diamond

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/myTrinotate.sqlite \
# --transcript_fasta /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BA/BA_bf.Trinity.fasta \
# --transdecoder_pep /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA/BA_bf.Trinity.fasta.transdecoder.pep \
# --trinotate_data_dir /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/trinotateDatDir \
# --run "swissprot_blastp swissprot_blastx pfam infernal EggnogMapper" --use_diamond

# #run TmHMM v2
# tmhmm --short /lustre04/scratch/janayfox/afFishRNA/trinotate/BN/BN_bf.Trinity.fasta.transdecoder.pep > tmhmm.v2.out
# tmhmm --short /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA/BA_bf.Trinity.fasta.transdecoder.pep > tmhmm.v2.out

#  #load into Trinotate 
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate/myTrinotate.sqlite \
# --LOAD_tmhmmv2 /lustre04/scratch/janayfox/afFishRNA/trinotate/tmhmm.v2.out 

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/myTrinotate.sqlite \
# --LOAD_tmhmmv2 /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/tmhmm.v2.out 


# #run signalp6
# signalp6 --fastafile /lustre04/scratch/janayfox/afFishRNA/trinotate/BN/BN_bf.Trinity.fasta.transdecoder.pep \
# --output_dir sigP6outdir_BN --format none --organism euk --mode fast

# signalp6 --fastafile /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/BA/BA_bf.Trinity.fasta.transdecoder.pep \
# --output_dir sigP6outdir_BA --format none --organism euk --mode fast

#load into Trinotate 
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate/myTrinotate.sqlite \
# --LOAD_signalp /lustre04/scratch/janayfox/afFishRNA/trinotate/sigP6outdir_BN/output.gff3

# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/myTrinotate.sqlite \
# --LOAD_signalp /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/sigP6outdir_BA/output.gff3

#generate report 
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate/myTrinotate.sqlite \
# --report --incl_pep --incl_trans > BN_trinotate.xls

singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
/usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate_BA/myTrinotate.sqlite \
--report --incl_pep --incl_trans > BA_trinotate.xls

