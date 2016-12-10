#!/usr/bin/env Rscript

# script collapses haplotypes (i.e. removes duplicate sequences)
# requires unaligned sequences in list, with no gap chars (otherwise gaps are considered differences)
# treats ambiguous bases (e.g. R, Y, N) as different 
# requires the 'stringr' and 'parallel' packages 
# returns the longest sequence if they are the otherwise the same, but different lengths

hapCollapse <- function(data, cores){#
    # sort by length
    data.ord <- data[order(mcmapply(length, data, SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=cores), decreasing=TRUE)]
    # make a copy to index later
    data.ord.copy <- data.ord
    # collapse into a strings
    data.ord <- mcmapply(FUN=function(x) paste(as.character(x), collapse=""), data.ord, SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=cores)
    # get the indices of the first match to each seq (the longest)
    ind <- unique(mcmapply(FUN=function(x) which(str_detect(string=data.ord, pattern=x) == TRUE)[1], data.ord, SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=cores))
    return(data.ord.copy[ind])
}#
