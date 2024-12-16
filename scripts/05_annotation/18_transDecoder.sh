#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Run Transdecoder on assemblies
### Author: Janay Fox
#######################################################

#extract long ORFs
./TransDecoder-v5.7.0/TransDecoder.LongOrfs -t ./readsBeforeRmoverrep/BN/BN_bf.Trinity.fasta
./TransDecoder-v5.7.0/TransDecoder.LongOrfs -t ./readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta

#predict coding regions
./TransDecoder-v5.7.0/TransDecoder.Predict -t ./readsBeforeRmoverrep/BN/BN_bf.Trinity.fasta
./TransDecoder-v5.7.0/TransDecoder.Predict -t ./readsBeforeRmoverrep/BA/BA_bf.Trinity.fasta
