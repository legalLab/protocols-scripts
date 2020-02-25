#!/usr/bin/env Rscript

# function to reverse complement DNA sequences in a dataframe
# Rupert Collins Feb 2020

# needs ape
library(ape)

# fun
revcomp_dna <- function(dnacol){
    dnas <- as.DNAbin(strsplit(dnacol,""))
    dnas.revcomp <- toupper(mapply(paste,collapse="",as.character(ape::complement(dnas)),SIMPLIFY=TRUE,USE.NAMES=FALSE))
    return(dnas.revcomp)
}

# examples
#ff <- tibble(acc_no=c("GB01","GB02","GB03"), species=c("human","cat","mouse"), sequence=c("TGGAAG","ATTTAC","ACCCAT"))
#ff %>% mutate(sequenceRevComp=revcomp_dna(dnacol=sequence))
