#!/usr/bin/env Rscript
# script collapses haplotypes - it requires a sequences in list, with no gap chars
# requires the parallel package
# returns the longest haplotype if they are the same but different lengths

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