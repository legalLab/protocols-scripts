#!/usr/bin/env Rscript

# this script removes all characters not a valid ACTG base from a DNAbin object. This includes:
# "N", "-", "?", "R", "Y", etc
# script gives a warning if these characters are inside the sequences, ie, are not alignment padding chars at the ends

#require("ape")
#require("stringr")

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



# another cleaning function
# slightly different - it can remove Ns from only the ends of the sequences, and you have to tell the function what you want to delete
# example:
# rm_missing(dna=cleantest, chars=c("?","-","n"), endOnlyNs=TRUE, print=TRUE)
# it's case insensitive too

rm_missing <- function(dna,chars,endOnlyNs,print){
# check if object is a DNAbin
if(class(dna)!="DNAbin"){stop("Please convert your data to a DNAbin object first")} else{
        dna.c <- as.character(as.list(dna))#convert dna to list and character
        dna.f <- sapply(dna.c, paste, collapse="")# cat dna
        chars[str_which(string=chars,pattern="\\?")] <- "\\?"# escape the question mark
            if(any(!str_detect(string=dna.f, pattern=regex("a|c|t|g", ignore_case=TRUE)))) {# finds sequences with no ACTGs
            warning("Some of the sequences (rows) have no ACTG nucleotides. I removed these and continued, therefore the length of the output will not be the same as the input.", call.=FALSE)
            dna.f <- dna.f[str_detect(string=dna.f, pattern=regex("a|c|t|g", ignore_case=TRUE))]# remove those and continue
            }
            if(endOnlyNs==TRUE) {
            chars[str_which(string=chars,pattern="n|N")] <- "%"# replace the n char with an % to remove the problem with empty vectors
            chars.l <- paste(chars, collapse="|")# create a list of chars
            chars.l <- paste(chars.l,"^n*|n*$",sep="|")
            } else {
                chars.l <- paste(chars, collapse="|")# create a list of chars
                }
        reg <- regex(pattern=chars.l, ignore_case=TRUE)# turn into a regex
        dna.r <- str_replace_all(string=dna.f, pattern=reg, replacement="")# replace the string
        names(dna.r) <- names(dna.f)# rename
        if(print==TRUE){print(as.list(dna.r))}
        dna.r <- as.DNAbin(str_split(string=dna.r, pattern=""))# convert back to a DNAbin
        names(dna.r) <- names(dna.f)# rename
        return(dna.r)
    }
}