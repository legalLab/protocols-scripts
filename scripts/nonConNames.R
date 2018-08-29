#!/usr/bin/env/ Rscript

# function 'nonConNames' to generate names of species that are the closest non-conspecifics to an individual
# Rupert Collins August 2018 
# see email to Pedro Senna

# fun
nonConNames <- function(distobj, sppVector, rmNA){
    distobj <- as.matrix(distobj)
    rownames(distobj) <- sppVector
    nonSpecNames <- list()
    for (i in 1:length(rownames(distobj))) {
        nonSpec <- rownames(distobj) != rownames(distobj)[i]
        nonSpecNames[[i]] <- paste(unique(names(which(distobj[nonSpec, i] == min(distobj[nonSpec, i],na.rm=rmNA), useNames=TRUE))), collapse="; ")
    }
    output <- unlist(nonSpecNames)
    output
}

# example
#library(ape)
#library(spider)
#data(anoteropsis)
#mat <- dist.dna(anoteropsis, model="raw", pairwise.deletion=TRUE)
#spp <- sapply(strsplit(dimnames(anoteropsis)[[1]], split="_"), function(x) paste(x[1], x[2], sep="_"))
#nonConNames(distobj=mat, sppVector=spp, rmNA=FALSE)
