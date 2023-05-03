#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-10:00
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

#BA
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db /lustre04/scratch/janayfox/afFishRNA/trinotate/myTrinotate.sqlite --init \
# --gene_trans_map /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BA/BA_bf.Trinity.fasta.gene_trans_map \
# --transcript_fasta /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BA/BA_bf.Trinity.fasta \
# --transdecoder_pep /lustre04/scratch/janayfox/afFishRNA/trinotate/BA/BA_bf.Trinity.fasta.transdecoder.pep

# #run sequence analyses
# singularity exec -e -B /lustre04/scratch/janayfox/afFishRNA trinotate.v4.0.0.simg \
# /usr/local/src/Trinotate/Trinotate --db myTrinotate.sqlite \
# --transcript_fasta /lustre04/scratch/janayfox/afFishRNA/trinity_output/bfRem_assembly/BN/BN_bf.Trinity.fasta \
# --transdecoder_pep /lustre04/scratch/janayfox/afFishRNA/trinotate/BN/BN_bf.Trinity.fasta.transdecoder.pep \
# --trinotate_data_dir /lustre04/scratch/janayfox/afFishRNA/trinotate/trinotateDatDir \
# --run "swissprot_blastp swissprot_blastx pfam infernal EggnogMapper" --use_diamond

#generate 

# #note I am running everything except signalp6 and tmhmmv2 which will need to be run separately 
# # singularity exec -e trinotate.v4.0.0.simg -B /lustre04/scratch/janayfox/afFishRNA --db trinotate_sqlite.db \
# # --transcript_fasta /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta \
# # --transdecoder_pep <insert> --trinotate_data_dir <insert> --run "swissprot_blastp swissprot_blastx pfam infernal EggnogMapper" \
# # --use_diamond


# #run TmHMM v2
tmhmm --short /lustre04/scratch/janayfox/afFishRNA/trinotate/BN/BN_bf.Trinity.fasta.transdecoder.pep > tmhmm.v2.out
 #load into Trinotate 
singularity exec -e trinotate.v4.0.0.simg -B /lustre04/scratch/janayfox/afFishRNA \
--db myTrinotate.sqlite --LOAD_tmhmmv2 tmhmm.v2.out #load into Trinotate 

# #run signalp6
# signalp6 --fastafile /lustre04/scratch/janayfox/afFishRNA/trinotate/BN/BN_bf.Trinity.fasta.transdecoder.pep --output_dir sigP6outdir --format none --organism euk --mode fast
# #load into Trinotate 
# singularity exec -e trinotate.v4.0.0.simg -B /lustre04/scratch/janayfox/afFishRNA \ 
# --db myTrinotate.sqlite --LOAD_signalp sigP6outdir/output.gff3

