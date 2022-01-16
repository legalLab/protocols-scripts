#!/usr/bin/env Rscript

# function to create a fasta file from tabular DNA sequences
# requires a dataframe and sequence data and names as columns
# Rupert Collins August 2018

# needs ape
#library(ape)

# fun
tab2fas <- function(df,seqcol,namecol){
    df <- as.data.frame(df)
    dtmp <- strsplit(as.character(df[[seqcol]]), "")
    names(dtmp) <- as.character(df[[namecol]])
    dat <- ape::as.DNAbin(dtmp)
    return(dat)
}

# examples
#ff <- data.frame(acc_no=c("GB01","GB02","GB03"), species=c("human","cat","mouse"), sequence=c("TGGAAG","ATTTAC","ACCCAT"))
#ff$names <- paste(ff$acc_no,ff$species,sep="_")
#gg <- tab2fas(df=ff, seqcol="sequence", namecol="names")
#write.dna(gg, file="gg.fas", format="fasta", colw=99999)


# FUNCTION TO CONVERT DNABIN OBJECT TO TABULAR
# works opposite to "tab2fas()"
# input can be aligned (matrix) or unaligned (list)
# strips out alignment characters
# requires "ape", "tidyr" and "stringr" packages
# returns a tibble
fas2tab <- function(dnas) {
    if(class(dnas) == "DNAbin") {
        dnas.list <- as.list(dnas)
        dnas.names <- names(dnas.list)
        dnas.char <- sapply(as.character(dnas.list),paste,collapse="")
        dnas.char.lower <- stringr::str_to_lower(dnas.char)
        dnas.clean <- stringr::str_replace_all(dnas.char.lower,"-|n|\\?","")
        dnas.tib <- tidyr::tibble(label=dnas.names,nucleotides=dnas.clean)
        return(dnas.tib)
    } else stop(writeLines("Object must be APE DNAbin format."))
}
