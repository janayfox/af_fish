######################################################################
### Goal: Extract results from Trinotate
### Author: Janay Fox
### R script
#####################################################################

## Set up ##
#read in data 
BN_blastx <- read.table("./data/trinotate/BN_uniprot_sprot.diamond.blastx.outfmt6")
BA_blastx <- read.table("./data/trinotate/BA_uniprot_sprot.diamond.blastx.outfmt6")

#add column names 
colnames(BN_blastx) <- c("q_seqid", "s_seqid", "p_ident", "length", "mismatch", 
                         "gap_open", "q_start", "q_end", "s_start", "s_end", 
                         "e_value", "bit_score")

colnames(BA_blastx) <- c("q_seqid", "s_seqid", "p_ident", "length", "mismatch", 
                         "gap_open", "q_start", "q_end", "s_start", "s_end", 
                         "e_value", "bit_score")
## Calculate stats ##
#identify unique proteins 
BN.prot <- unique(BN_blastx$s_seqid)
BA.prot <- unique(BA_blastx$s_seqid)

#filter to only matches that are 80% 
BN_blastx_80 <- BN_blastx[BN_blastx$p_ident >= 80,]
BA_blastx_80 <- BA_blastx[BA_blastx$p_ident >= 80,]

#identify unique proteins 
BN.prot_80 <- unique(BN_blastx_80$s_seqid)
BA.prot_80 <- unique(BA_blastx_80$s_seqid)
