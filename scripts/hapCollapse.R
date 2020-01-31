#!/usr/bin/env Rscript
library("ape")
library("stringr")

# script collapses haplotypes (i.e. removes duplicate sequences)
# if `collapseSubstrings==TRUE` the function will consider shorter but identical sequences as the same haplotype and collapse them, returning the longest sequence
# if `collapseSubstrings==FALSE` the function will consider shorter but identical sequences as different haplotypes and will keep them.
# if `clean==TRUE` the function will remove all non "A,C,T,G" positions such as "N,?,-" but also ambiguity codes such as "R,Y".
# if `clean==FALSE` the function will treat the data as is, and will not remove any bases
# requires the 'stringr' and 'ape' packages 
# 

hapCollapse <- function(data,collapseSubstrings,clean){
    source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/clean_dna.R")
    if(clean==TRUE){
    data <- clean_dna(data)
    }
    else{
    data <- as.list(data)
    }
    if(collapseSubstrings==TRUE){
    # sort by length
    data.ord <- data[order(mapply(length, data, SIMPLIFY=TRUE, USE.NAMES=FALSE), decreasing=TRUE)]
    # make a copy to index later
    data.ord.copy <- data.ord
    # collapse into a strings
    data.ord <- mapply(FUN=function(x) paste(x,collapse=""), as.character(data.ord), SIMPLIFY=TRUE, USE.NAMES=FALSE)
    # get the indices of the first match to each seq (the longest)
    ind <- unique(mapply(FUN=function(x) which(str_detect(string=data.ord, pattern=x) == TRUE)[1], data.ord, SIMPLIFY=TRUE, USE.NAMES=FALSE))
    return(data.ord.copy[ind])
    }
    else{
    data.copy <- data
    data.char <- mapply(FUN=function(x) paste(x,collapse=""), as.character(data), SIMPLIFY=TRUE, USE.NAMES=FALSE)
    dups <- duplicated(data)
    data.keep <- data.copy[!dups]
    return(data.keep)
    }
}
