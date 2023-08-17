# author: dlueckin
# date: Mon Jan 16 11:38:12 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# import blast result -----------------------------------------------------

blast_out <- fread("top20/gtas/blastp_result.out")
names(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast_out <- blast_out %>% 
    filter(pident > 50 & evalue < 10^-5)


# read proteins, and select relevant ones ---------------------------------

all_proteins <- read.fasta(file = "top20/gtas/combined_top20_proteins.faa", seqtype = "AA")
relevant_proteins <- all_proteins[which(names(all_proteins) %in% blast_out$sseqid)]



# write out relevant_proteins in batches of 30 for interpro ---------------

n = 1
for(i in seq(1, length(relevant_proteins), by = 30)){
    write.fasta(relevant_proteins[i:(i+29)],
                names = names(relevant_proteins)[i:(i+29)],
                file.out = paste0("top20/gtas/interpro/input_proteins/input_proteins_", n, ".faa"))
    n = n+1
}













