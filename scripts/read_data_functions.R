#!/usr/bin/env Rscript

# functions to read various genetic data formats formats
# depends on adegenet to convert df objects (read_genepop)


################################
#' @title read_genepop
#' @description read .gen formated datasets
#' @description adapted read.genepop function (adegenet package) to allow reading population information
#' @author Tomas Hrbek February 2021
#'
#' @param file -> gen formatted file
#' @param ncode -> length of word coding each allele (integer)
#' @param quiet -> displays runtime information
#' @export nothing
#' @return gen class object
#'
#' @details
#' This function reads a gen formatted infile and returns it as gen class object
#' Unlike read.genepop it stores population information and returns it
#' Unlike read.genepop it does not throw an exception if infile does not end with .gen
#' 
#' @example
#' read_genepop(file = my_infile, ncode = 2L, quiet = FALSE)
#' read_genepop(my_infile, node = 2L)
#' read_genepop(my_infile)
#'

read_genepop <- function(file, ncode = 2L, quiet = FALSE) 
{
#  if (!require(adegenet)) {
#    install.packages("adegenet")
#  }
# permit files with an extension other than .gen
#  if (toupper(.readExt(file)) != "GEN") 
#    stop("File extension .gen expected")
  if (!quiet) 
    cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")
  prevcall <- match.call()
  txt <- scan(file, sep = "\n", what = "character", quiet = TRUE)
  if (!quiet) 
    cat("\nFile description: ", txt[1], "\n")
  txt <- txt[-1]
  txt <- gsub("\t", " ", txt)
  NA.char <- paste(rep("0", ncode), collapse = "")
  locinfo.idx <- 1:(min(grep("POP", toupper(txt))) - 1)
  locinfo <- txt[locinfo.idx]
  locinfo <- paste(locinfo, collapse = ",")
  loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
  loc.names <- trimws(loc.names)
  nloc <- length(loc.names)
  txt <- txt[-locinfo.idx]
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)", toupper(txt))
  npop <- length(pop.idx)
  nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
  splited <- nocomma[which(!nocomma %in% pop.idx)]
  if (length(splited) > 0) {
    for (i in sort(splited, decreasing = TRUE)) {
      txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
    }
    txt <- txt[-splited]
  }
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)", toupper(txt))
  txt[length(txt) + 1] <- "POP"
  nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)", 
                          toupper(txt))) - 1
  pop <- factor(rep(1:npop, nind.bypop))
  temp <- sapply(1:length(txt[pop.idx]), function(i) strsplit(txt[pop.idx][i], "[[:space:]]+"))
  pop.names <- vapply(temp, function(x) x[2], character(1))
  pop.names <- trimws(pop.names)
# keep names of populations
# pop.names <- replace(pop.names, is.na(pop.names),"pop")
  levels(pop) <- pop.names
  txt <- txt[-c(pop.idx, length(txt))]
  temp <- sapply(1:length(txt), function(i) strsplit(txt[i], 
                                                     ","))
  ind.names <- vapply(temp, function(x) x[1], character(1))
  ind.names <- trimws(ind.names)
  vec.genot <- vapply(temp, function(x) x[2], character(1))
  vec.genot <- trimws(vec.genot)
  X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), 
              ncol = nloc, byrow = TRUE)
  if (any(duplicated(ind.names))) {
    rownames(X) <- .genlab("", nrow(X))
    warning("Duplicate individual names detected. Coercing them to be unique.")
  }
  else {
    rownames(X) <- ind.names
  }
  colnames(X) <- loc.names
  if (!all(unique(nchar(X)) == (ncode * 2))) 
    stop(paste("some alleles are not encoded with", ncode, 
               "characters\nCheck 'ncode' argument"))
  res <- adegenet::df2genind(X = X, pop = as.character(pop), ploidy = 2, 
                   ncode = ncode, NA.char = NA.char)
  res@call <- prevcall
  if (!quiet) 
    cat("\n...done.\n\n")
  return(res)
}