#!/usr/bin/env Rscript

# improvements to existing adegenet functions
# read genepop format with population information and with other extensions (read_genepop)
# run PCA on genlight data and output prcomp object for downstream analysis and plotting (snp_pca)


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


################################
#' @title snp_pca
#' @description run PCA on SNP data in genlight format
#' @description adapted glPca function (adegenet package) to output prcomp class object which includes all prcomp class object fields
#' @author Tomas Hrbek February 2021
#'
#' @param file -> genlight formatted file
#' @param center -> subtract variable mean; only effect of variance remains (logical)
#' @param scale -> z-score normalization to mean = 0, sd = 1 (logical)
#' @param nf -> number of factors (PCs) to calculate (integer)
#' @param loadings -> retain rotation (loci x PCs) (logical)
#' @param other parameters from glPca
#' @export nothing
#' @return prcomp class object
#'
#' @details
#' This function runs a PCA on genlight formated data and returns it as prcomp class object
#' Unlike glPca returns prcomp object with all information found in prcomp objects
#' Unlike glPca result can be plotted using autoplot() and ggplot() functions
#' Unlike glPca scales data by default scale = TRUE
#' 
#' @example
#' snp_pca(file = my_genlight, center = TRUE, scale = FALSE, nf = NULL)
#' snp_pca(my_genlight, nf = 5)
#' snp_pca(my_genlight)
#'

snp_pca <- function (file, center = TRUE, scale = TRUE, nf = NULL, loadings = TRUE, 
                     alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE, n.cores = NULL, 
                     returnDotProd = FALSE, matDotProd = NULL) {
  if (!inherits(file, "genlight")) 
    stop("file is not a genlight object")
  if (center) {
    vecMeans <- glMean(file, alleleAsUnit = alleleAsUnit)
    if (any(is.na(vecMeans))) 
      stop("NAs detected in the vector of means")
  }
  if (scale) {
    vecVar <- glVar(file, alleleAsUnit = alleleAsUnit)
    if (any(is.na(vecVar))) 
      stop("NAs detected in the vector of variances")
  }
  myPloidy <- ploidy(file)
  if (is.null(matDotProd)) {
    if (!useC) {
      if (parallel && is.null(n.cores)) {
        n.cores <- parallel::detectCores()
      }
      if (!center & !scale) {
        dotProd <- function(a, b, ploid.a, ploid.b) {
          a <- as.integer(a)/ploid.a
          a[is.na(a)] <- 0
          b <- as.integer(b)/ploid.b
          b[is.na(b)] <- 0
          return(sum(a * b, na.rm = TRUE))
        }
      }
      if (center & !scale) {
        dotProd <- function(a, b, ploid.a, ploid.b) {
          a <- as.integer(a)/ploid.a
          a[is.na(a)] <- vecMeans[is.na(a)]
          b <- as.integer(b)/ploid.b
          b[is.na(b)] <- vecMeans[is.na(b)]
          return(sum((a - vecMeans) * (b - vecMeans), na.rm = TRUE))
        }
      }
      if (!center & scale) {
        dotProd <- function(a, b, ploid.a, ploid.b) {
          a <- as.integer(a)/ploid.a
          a[is.na(a)] <- 0
          b <- as.integer(b)/ploid.b
          b[is.na(b)] <- 0
          return(sum((a * b)/vecVar, na.rm = TRUE))
        }
      }
      if (center & scale) {
        dotProd <- function(a, b, ploid.a, ploid.b) {
          a <- as.integer(a)/ploid.a
          a[is.na(a)] <- vecMeans[is.na(a)]
          b <- as.integer(b)/ploid.b
          b[is.na(b)] <- vecMeans[is.na(b)]
          return(sum(((a - vecMeans) * (b - vecMeans))/vecVar, na.rm = TRUE))
        }
      }
      allComb <- combn(1:nInd(file), 2)
      if (parallel) {
        allProd <- unlist(parallel::mclapply(1:ncol(allComb), function(i) dotProd(file@gen[[allComb[1, i]]], file@gen[[allComb[2, i]]], myPloidy[allComb[1, i]], myPloidy[allComb[2, i]]), mc.cores = n.cores, mc.silent = TRUE, mc.cleanup = TRUE, mc.preschedule = FALSE))
      }
      else {
        allProd <- unlist(lapply(1:ncol(allComb), function(i) dotProd(file@gen[[allComb[1, i]]], file@gen[[allComb[2, i]]], myPloidy[allComb[1, i]], myPloidy[allComb[2, i]])))
      }
      allProd <- allProd/nInd(file)
      attr(allProd, "Size") <- nInd(file)
      attr(allProd, "Diag") <- FALSE
      attr(allProd, "Upper") <- FALSE
      class(allProd) <- "dist"
      allProd <- as.matrix(allProd)
      if (parallel) {
        temp <- unlist(parallel::mclapply(1:nInd(file), function(i) dotProd(file@gen[[i]], file@gen[[i]], myPloidy[i], myPloidy[i]), mc.cores = n.cores, mc.silent = TRUE, mc.cleanup = TRUE, mc.preschedule = FALSE))/nInd(file)
      }
      else {
        temp <- unlist(lapply(1:nInd(file), function(i) dotProd(file@gen[[i]], file@gen[[i]], myPloidy[i], myPloidy[i])))/nInd(file)
      }
      diag(allProd) <- temp
    }
    else {
      allProd <- glDotProd(file, center = center, scale = scale, alleleAsUnit = alleleAsUnit, parallel = parallel, n.cores = n.cores)/nInd(file)
    }
  }
  else {
    if (!all(dim(matDotProd) == nInd(file))) 
      stop("matDotProd has wrong dimensions.")
    allProd <- matDotProd
  }
  eigRes <- eigen(allProd, symmetric = TRUE, only.values = FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop = FALSE]
  if (is.null(nf)) {
    barplot(eigRes$values, main = "Eigenvalues", col = heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(con = getOption("adegenet.testcon"), 
                               n = 1))
  }
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(eigRes$values > 1e-10))
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(file))
  res$x <- sweep(eigRes$vectors[, 1:nf, drop = FALSE], 2, sqrt(eigRes$values[1:nf]), FUN = "*")
  if (loadings) {
    if (scale) {
      vecSd <- sqrt(vecVar)
    }
    res$rotation <- matrix(0, nrow = nLoc(file), ncol = nf)
    for (k in 1:nInd(file)) {
      temp <- as.integer(file@gen[[k]])/myPloidy[k]
      if (center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      }
      else {
        temp[is.na(temp)] <- 0
      }
      if (scale) {
        temp <- temp/vecSd
      }
      res$rotation <- res$rotation + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop = FALSE]
    }
    res$rotation <- res$rotation/nInd(file)
    res$rotation <- sweep(res$rotation, 2, sqrt(eigRes$values[1:nf]), 
                          FUN = "/")
    #res$rotation <- res$rotation[rowSums(res$rotation) != 0,]
  }
  colnames(res$x) <- paste("PC", 1:nf, sep = "")
  if (!is.null(indNames(file))) {
    rownames(res$x) <- indNames(file)
  }
  else {
    rownames(res$x) <- 1:nInd(file)
  }
  if (!is.null(res$rotation)) {
    colnames(res$rotation) <- paste("PC", 1:nf, sep = "")
    if (!is.null(locNames(file)) & !is.null(alleles(file))) {
      rownames(res$rotation) <- paste(locNames(file), alleles(file), sep = ".")
    }
    else {
      rownames(res$rotation) <- 1:nLoc(file)
    }
  }
  if (returnDotProd) {
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(file)
  }
  #res$rotation <- as(res$rotation, "matrix")
  res$sdev <- apply(res$x, 2, sd)
  res$center <- glMean(file, alleleAsUnit = alleleAsUnit)
  res$scale <- sqrt(glVar(file, alleleAsUnit = alleleAsUnit))
  res$pop <- pop(file)
  res$call <- match.call()
  class(res) <- "prcomp"
  return(res)
}

