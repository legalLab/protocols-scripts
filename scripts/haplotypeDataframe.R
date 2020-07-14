#!/usr/bin/env Rscript

# load libs and functions
library("tidyverse")
library("ape")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/hapCollapse.R")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/clean_dna.R")

# load dna
dat.all <- clean_dna(read.FASTA("https://raw.githubusercontent.com/legalLab/publications/master/Bittencourt_et_al_2019/P_trigonatus_250718_ali.fas"))

# collapse to haplotypes
dat.haps <- hapCollapse(clean_dna(dna=dat.all),cores=2)

# convert to all to character
dat.haps.char <- lapply(dat.haps, function(x) paste(x, collapse=""))
dat.all.char <- lapply(dat.all, function(x) paste(x, collapse=""))
dat.daughters.char <- dat.all.char[which(!names(dat.all.char) %in% names(dat.haps.char))]

# detect strings for all daughters of each unique haplotype
seqs.in <- mapply(FUN=function(x) which(str_detect(string=dat.haps.char, pattern=x)==TRUE), dat.daughters.char, SIMPLIFY=TRUE, USE.NAMES=FALSE)
seqs.in <- mapply(function(x) names(dat.haps.char[x]), seqs.in)
names(seqs.in) <- names(dat.daughters.char)

# turn into dataframe
pair.list <- unnest(enframe(seqs.in,name="daughter",value="mother"))

# view
print(pair.list)
