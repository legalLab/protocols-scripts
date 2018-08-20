#!/usr/bin/env Rscript

# function to create a fasta file from tabular DNA sequences
# requires a dataframe and sequence data and names as columns
# Rupert Collins August 2018

# needs ape
library(ape)

# fun
tab2fas <- function(df,seqcol,namecol){
    dtmp <- strsplit(as.character(df[[seqcol]]), "")
    names(dtmp) <- as.character(df[[namecol]])
    dat <- as.DNAbin(dtmp)
    return(dat)
}

# examples
#ff <- data.frame(acc_no=c("GB01","GB02","GB03"), species=c("human","cat","mouse"), sequence=c("TGGAAG","ATTTAC","ACCCAT"))
#ff$names <- paste(ff$acc_no,ff$species,sep="_")
#gg <- tab2fas(df=ff, seqcol="sequence", namecol="names")
#write.dna(gg, file="gg.fas", format="fasta", colw=99999)
