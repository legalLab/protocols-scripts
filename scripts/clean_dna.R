#!/usr/bin/env Rscript

# this script removes all characters not a valid ACTG base from a DNAbin object. This includes:
# "N", "-", "?", "R", "Y", etc
# script gives a warning if these characters are inside the sequences, ie, are not alignment padding chars at the ends

require("ape")
require("stringr")

clean_dna <- function(dna){#
    # convert to a list
    dat <- as.list(dna)
    # convert to character
    datc <- as.character(dat)
    # make a warning
        # collapse text into vector
        cdat <- lapply(datc, function(x) paste(x, collapse=""))
        # find all internal missing data
        res <- lapply(cdat, function(x) str_detect(string=x, pattern="[actg][^actg]+[actg]"))
        # get names
        errs <- names(which(res!=FALSE))
        # if else warning
        if(length(errs >= 1)){print(errs); warning("You have missing data ('N','-' '?') or ambiguity inside your sequence, i.e. not padding the ends, and this may have unintended consequences later, as they have now been removed! The names of the samples are above.")}
    # get positions of all non bases
    inds <- lapply(datc, function(x) grep("a|c|t|g", x, ignore.case=TRUE))
    # match positions and remove  
    dlimp <- mapply(function(x,y) x[y], x=datc, y=inds, SIMPLIFY=FALSE)
    # convert to dnabin
    dbin <- as.DNAbin(dlimp)
    return(as.list(dbin))
}#