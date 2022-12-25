#!/usr/bin/env Rscript

# functions to convert vcf data to other formats
# depends on vcfR (should be loaded since vcfR object [vcf file] is passed to these functions)
# depends on dplyr (many functions)
# depends on ape (vcf2nexus, vcf2snapp)
# depends on glue (vcf2fineRadStructure, vcf2treemix)
# depends on adegenet (vcf2genlight)


################################
#' @title vcf2migrate
#' @description converts vcfR format data to MigrateN infile
#' @description adapted vcfR2migrate function (vcfR package) to allow for inclusion of missing data in Migrate output
#' @author Tomas Hrbek August 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in MigrateN infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (MigrateN infile)
#' @param method -> classic or het format
#' @export MigrateN infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a MigrateN formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2migrate(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "MigrateN_infile.txt", method = "N")
#' vcf2migrate(my_vcf, ind_pop, keepers, out_file = "MigrateN_infile.txt")
#' vcf2migrate(my_vcf, ind_pop, keepers)
#'

vcf2migrate <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE, 
                         out_file = "migrateN_infile.txt", method = "N")
{
  method <- match.arg(method, c("N", "H"), several.ok = FALSE)
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  
  if (method == "N") {
    myHeader <- c("N", length(vcf_list), nrow(vcf_list[[1]]))
    pop_list <- vector(mode = "list", length = length(vcf_list))
    names(pop_list) <- names(vcf_list)
    for (i in 1:length(vcf_list)) {
      gt <- extract.gt(vcf_list[[i]], return.alleles = T, convertNA = T) #convertNA not working here
      gt[gt == "."] <- "?/?"
      allele1 <- apply(gt, MARGIN = 2, function(x) {
        substr(x, 1, 1)
      })
      rownames(allele1) <- NULL
      allele1 <- t(allele1)
      rownames(allele1) <- paste(rownames(allele1), "_1", 
                                 sep = "")
      allele2 <- apply(gt, MARGIN = 2, function(x) {
        substr(x, 3, 3)
      })
      rownames(allele2) <- NULL
      allele2 <- t(allele2)
      rownames(allele2) <- paste(rownames(allele2), "_2", 
                                 sep = "")
      pop_list[[i]][[1]] <- allele1
      pop_list[[i]][[2]] <- allele2
    }
    write(myHeader, file = out_file, ncolumns = length(myHeader), 
          sep = "\t")
    write(rep(1, times = ncol(pop_list[[1]][[1]])), file = out_file, 
          ncolumns = ncol(pop_list[[1]][[1]]), append = TRUE, 
          sep = "\t")
    for (i in 1:length(pop_list)) {
      popName <- c(2 * nrow(pop_list[[i]][[1]]), names(pop_list)[i])
      write(popName, file = out_file, ncolumns = length(popName), 
            append = TRUE, sep = "\t")
      for (j in 1:ncol(pop_list[[i]][[1]])) {
        utils::write.table(pop_list[[i]][[1]][, j], file = out_file, 
                           append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, 
                           col.names = FALSE)
        utils::write.table(pop_list[[i]][[2]][, j], file = out_file, 
                           append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, 
                           col.names = FALSE)
      }
    }
  }
  else if (method == "H") {
    myHeader <- c("H", length(vcf_list), nrow(vcf_list[[1]]))
    pop_list <- vector(mode = "list", length = length(vcf_list))
    names(pop_list) <- names(vcf_list)
    for (i in 1:length(vcf_list)) {
      myMat <- matrix(nrow = nrow(vcf_list[[i]]), ncol = 6)
      var_info <- as.data.frame(vcf_list[[i]]@fix[, 1:2, 
                                                  drop = FALSE])
      var_info$mask <- TRUE
      gt <- extract.gt(vcf_list[[i]])
      popSum <- .gt_to_popsum(var_info = var_info, gt = gt)
      myMat[, 1] <- paste(vcf_list[[i]]@fix[, "CHROM"], 
                          vcf_list[[i]]@fix[, "POS"], sep = "_")
      myMat[, 2] <- vcf_list[[i]]@fix[, "REF"]
      myMat[, 4] <- vcf_list[[i]]@fix[, "ALT"]
      myMat[, 3] <- unlist(lapply(strsplit(as.character(popSum$Allele_counts), 
                                           split = ",", fixed = TRUE), function(x) {
                                             x[1]
                                           }))
      myMat[, 3][is.na(myMat[, 3])] <- 0
      myMat[, 5] <- unlist(lapply(strsplit(as.character(popSum$Allele_counts), 
                                           split = ",", fixed = TRUE), function(x) {
                                             x[2]
                                           }))
      myMat[, 5][is.na(myMat[, 5])] <- 0
      myMat[, 6] <- as.numeric(myMat[, 3]) + as.numeric(myMat[, 
                                                              5])
      pop_list[[i]] <- myMat
    }
    write(myHeader, file = out_file, ncolumns = length(myHeader), 
          sep = "\t")
    for (i in 1:length(pop_list)) {
      popName <- c(pop_list[[i]][1, 6], names(pop_list[i]))
      write(popName, file = out_file, ncolumns = length(popName), 
            append = TRUE, sep = "\t")
      utils::write.table(pop_list[[i]], file = out_file, 
                         append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, 
                         col.names = FALSE)
    }
  }
  else {
    stop("You should never get here!")
  }
  return(invisible(NULL))
}


################################
#' @title vcf2arlequin
#' @description converts vcfR format data to Arlequin infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Arlequin infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Arlequin infile)
#' @export Arlequin infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Arlequin formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2arlequin(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Arlequin_infile.arp")
#' vcf2arlequin(my_vcf, ind_pop, keepers, out_file = "Arlequin_infile.arp")
#' vcf2arlequin(my_vcf, ind_pop, keepers)
#'

vcf2arlequin <-function (vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "arlequin.arp") 
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  pop_list <- vector(mode = "list", length = length(vcf_list))
  names(pop_list) <- names(vcf_list)
  
  for (i in 1:length(vcf_list)) {
    gt <- extract.gt(vcf_list[[i]], return.alleles = T, convertNA = T) #convertNA not working here
    gt[gt == "."] <- "?/?"
    allele1 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 1, 1)
    })
    rownames(allele1) <- NULL
    allele1 <- t(allele1)
    rownames(allele1) <- paste(rownames(allele1), "_1", 
                               sep = "")
    allele2 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 3, 3)
    })
    rownames(allele2) <- NULL
    allele2 <- t(allele2)
    rownames(allele2) <- paste(rownames(allele2), "_2", 
                               sep = "")
    pop_list[[i]][[1]] <- allele1
    pop_list[[i]][[2]] <- allele2
  }
  
  write("[Profile]", file = out_file)
  write("", file = out_file, append = TRUE)
  write("Title = 'Generated by vcf2arlequin.R'", file = out_file, append = TRUE)
  write(paste("NbSamples = ", length(vcf_list), sep = ""), file = out_file, append = TRUE)
  write("GenotypicData = 1", file = out_file, append = TRUE)
  write("LocusSeparator = WHITESPACE", file = out_file, append = TRUE)
  write("GameticPhase = 0", file = out_file, append = TRUE)
  write("MissingData = '?'", file = out_file, append = TRUE)
  write("DataType = STANDARD", file = out_file, append = TRUE)
  write("", file = out_file, append = TRUE)
  write("[Data]", file = out_file, append = TRUE)
  write("[[Samples]]", file = out_file, append = TRUE)
  write("", file = out_file, append = TRUE)
  write(paste("#There are ", nrow(vcf), " SNPs", sep = ""), file = out_file, append = TRUE)
  write("", file = out_file, append = TRUE)
  
  for (i in 1:length(pop_list)) {
    write(paste("SampleName = ", "'", names(pop_list)[i], "'", sep = ""), file = out_file, append = TRUE)
    write(paste("SampleSize = ", nrow(pop_list[[i]][[1]]), sep = ""), file = out_file, append = TRUE)
    write("SampleData={", file = out_file, append = TRUE)
    
    for (j in 1:nrow(pop_list[[i]][[1]])) {
      utils::write.table(t(c(names(pop_list[[i]][[1]][j, 1]), "\t1\t", pop_list[[i]][[1]][j, ])), file = out_file, 
                         append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, 
                         col.names = FALSE)
      utils::write.table(t(c("\t\t\t\t\t", pop_list[[i]][[2]][j, ])), file = out_file, 
                         append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, 
                         col.names = FALSE)
    }
    write("}", file = out_file, append = TRUE)
    write("", file = out_file, append = TRUE)
  }
  
  write("[[Structure]]", file = out_file, append = TRUE)
  write("StructureName = 'One Group'", file = out_file, append = TRUE)
  write("NbGroups = 1", file = out_file, append = TRUE)
  write("", file = out_file, append = TRUE)
  write("Group = {",  file = out_file, append = TRUE)
  for (i in 1:length(names(pop_list))) {
    write(paste("\t\t", "\"", names(pop_list)[i], "\"", sep = ""),  file = out_file, append = TRUE)
  }
  write("}",  file = out_file, append = TRUE)
  
  return(invisible(NULL))
}


################################
#' @title vcf2structure
#' @description converts vcfR format data to Structure or FastStructure infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Structure infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Structure infile)
#' @param method -> Structure or FastStructure format
#' @export Structure infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Structure or FastStructure formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2structure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Structure_infile.str", method = "S")
#' vcf2structure(my_vcf, ind_pop, keepers, out_file = "Structure_infile.str")
#' vcf2structure(my_vcf, ind_pop, keepers)
#'

vcf2structure <-function (vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "structure.str", method = "S") 
{
  method <- match.arg(method, c("S", "F"), several.ok = FALSE)
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  pop_list <- vector(mode = "list", length = length(vcf_list))
  names(pop_list) <- names(vcf_list)
  
  for (i in 1:length(vcf_list)) {
    gt <- extract.gt(vcf_list[[i]], return.alleles = F, convertNA = T) #convertNA not working here
    gt[is.na(gt)] <- "?/?"
    allele1 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 1, 1)
    })
    rownames(allele1) <- NULL
    allele1 <- t(allele1)
    allele1[allele1 == "?"] <- "-9"
    rownames(allele1) <- paste(rownames(allele1), "_1", 
                               sep = "")
    allele2 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 3, 3)
    })
    rownames(allele2) <- NULL
    allele2 <- t(allele2)
    allele2[allele2 == "?"] <- "-9"
    rownames(allele2) <- paste(rownames(allele2), "_2", 
                               sep = "")
    pop_list[[i]][[1]] <- allele1
    pop_list[[i]][[2]] <- allele2
  }
  
  if (file.exists(out_file)) {
    file.remove(out_file)
  }
  
  # default output Structure, alternate output FastStructure
  if (method == "S") {
    for (i in 1:length(pop_list)) {
      for (j in 1:nrow(pop_list[[i]][[1]])) {
        utils::write.table(t(c(names(pop_list[[i]][[1]][j, 1]), i, pop_list[[i]][[1]][j, ])), file = out_file, 
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, 
                           col.names = FALSE)
        utils::write.table(t(c(names(pop_list[[i]][[2]][j, 1]), i, pop_list[[i]][[2]][j, ])), file = out_file, 
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, 
                           col.names = FALSE)
      }
    }
  } else if (method == "F") {
    fill <- rep(c(0), 4)
    for (i in 1:length(pop_list)) {
      for (j in 1:nrow(pop_list[[i]][[1]])) {
        utils::write.table(t(c(names(pop_list[[i]][[1]][j, 1]), i, fill, pop_list[[i]][[1]][j, ])), file = out_file, 
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, 
                           col.names = FALSE)
        utils::write.table(t(c(names(pop_list[[i]][[2]][j, 1]), i, fill, pop_list[[i]][[2]][j, ])), file = out_file, 
                           append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, 
                           col.names = FALSE)
      }
    }
  }
  
  return(invisible(NULL))
}


################################
#' @title vcf2genepop
#' @description converts vcfR format data to Genepop infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek January 2021
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Genepop infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Genepop infile)
#' @export Genepop infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Genepop formatted input file
#' This function labels populations. To read labeled populations use "read_genepop" function
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#' 01, 02, 03, 04 is 'A', 'C', 'G', 'T'
#'
#' @example
#' vcf2genepop(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Genepop_infile.gen")
#' vcf2genepop(my_vcf, ind_pop, keepers, out_file = "Genepop_infile.gen")
#' vcf2genepop(my_vcf, ind_pop, keepers)
#'

vcf2genepop <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE, 
                         out_file = "genepop_infile.txt")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  pop_list <- vector(mode = "list", length = length(vcf_list))
  names(pop_list) <- names(vcf_list)
  
  gt <- extract.gt(vcf, return.alleles = F, convertNA = T)
  gt <- t(gt)
  
  write("Title = 'Generated by vcf2genpop.R'", file = out_file)
  suppressWarnings(write.table(gt[0,], file = out_file, quote = FALSE, sep = ", ", col.names = TRUE, append = TRUE))
  
  for (i in 1:length(vcf_list)) {
    gt <- extract.gt(vcf_list[[i]], return.alleles = T, convertNA = T) #convertNA not working here
    gt[gt == "."] <- "0000"
    gt[gt == "A/A" | gt == "A|A"] <- "0101"
    gt[gt == "A/C" | gt == "C/A" | gt == "A|C" | gt == "C|A"] <- "0102"
    gt[gt == "A/G" | gt == "G/A" | gt == "A|G" | gt == "G|A"] <- "0103"
    gt[gt == "A/T" | gt == "T/A" | gt == "A|T" | gt == "T|A"] <- "0104"
    gt[gt == "C/C" | gt == "C|C"] <- "0202"
    gt[gt == "C/G" | gt == "G/C" | gt == "C|G" | gt == "G|C"] <- "0203"
    gt[gt == "C/T" | gt == "T/C" | gt == "C|T" | gt == "T|C"] <- "0204"
    gt[gt == "G/G" | gt == "G|G"] <- "0303"
    gt[gt == "G/T" | gt == "T/G" | gt == "G|T" | gt == "T|G"] <- "0304"
    gt[gt == "T/T" | gt == "T|T"] <- "0404"
    gt <- t(gt)
    write(paste("pop ", names(pop_list)[i], sep = ""), file = out_file, append = TRUE)
    utils::write.table(cbind(sep = ',', gt), file = out_file, quote = FALSE, sep = " ", col.names = FALSE, append = TRUE)
  }
  
  return(invisible(NULL))
}


################################
#' @title vcf2genlight
#' @description converts vcfR format data to Genlight infile
#' @description in part based on vcfR2genlight function (vcfR package)
#' @author Tomas Hrbek December 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Genlight infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @export Genlight object
#' @return Genlight object
#'
#' @details
#' This function converts the vcfR object to a Genlight formatted input file
#' This function labels populations.
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2genlight(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, ploidy = 2, inc_missing = TRUE, out_file = "Genlight_infile.txt")
#' vcf2genlight(my_vcf, ind_pop, keepers, out_file = "Genlight_infile.txt")
#' vcf2genlight(my_vcf, ind_pop, keepers)

vcf2genlight <- function (vcf, ind_pop, keep_pop, ploidy = 2, inc_missing = TRUE, 
                          out_file = "genlight_infile.R")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_sub_pops(vcf, ind_pop, keep_pop)
  if (any(is.na(getID(vcf2)))) {
    vcf2 <- addID(vcf2)
  }
  gt <- extract.gt(vcf2, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "NA"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "2"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "1"
  
  # create genelight object
  suppressWarnings(x <- as(t(gt), "genlight"))
  
  adegenet::chromosome(x) <- getCHROM(vcf2)
  adegenet::position(x) <- getPOS(vcf2)
  adegenet::pop(x) <- ind_pop[ind_pop %in% keep_pop]
  adegenet::ploidy(x) <- ploidy
  
  save(x, file = out_file)
  
  return(x)
}


################################
#' @title vcf2bayescan
#' @description converts vcfR format data to Bayescan infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek November 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Bayescan infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Bayescan infile)
#' @export Bayescan infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Bayescan formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2bayescan(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Bayescan_infile.txt")
#' vcf2bayescan(my_vcf, ind_pop, keepers, out_file = "Bayescan_infile.txt")
#' vcf2bayescan(my_vcf, ind_pop, keepers)
#'

vcf2bayescan <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE, 
                         out_file = "bayescan_infile.txt")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  pop_list <- vector(mode = "list", length = length(vcf_list))
  names(pop_list) <- names(vcf_list)
  
  for (i in 1:length(vcf_list)) {
    gt <- extract.gt(vcf_list[[i]], return.alleles = F, convertNA = T) #convertNA not working here
    allele1 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 1, 1)
    })
    rownames(allele1) <- c(1:nrow(allele1))
    allele2 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 3, 3)
    })
    rownames(allele2) <- c(1:nrow(allele1))
    
    pop_list[[i]][[3]] <- apply(cbind(allele1, allele2), 1, function(x){sum(x == 0, na.rm = TRUE)})
    pop_list[[i]][[4]] <- apply(cbind(allele1, allele2), 1, function(x){sum(x == 1, na.rm = TRUE)})
    pop_list[[i]][[1]] <- pop_list[[i]][[3]] + pop_list[[i]][[4]]
    pop_list[[i]][[2]] <- rep(2, length(pop_list[[i]][[1]]))
  }
  
  write(paste("[loci]=", nrow(vcf), sep = ""), file = out_file)
  write("", file = out_file, append = TRUE)
  write(paste("[populations]=", length(pop_list), sep = ""), file = out_file, append = TRUE)
  
  for (i in 1:length(pop_list)) {
    write("", file = out_file, append = TRUE)
    write(paste("[pop]=", i, sep = ""), file = out_file, append = TRUE)
    utils::write.table(pop_list[[i]], file = out_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE, append = TRUE)
  }
  
  return(invisible(NULL))
}


################################
#' @title vcf2treemix
#' @description converts vcfR format data to TreeMix infile
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek November 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in TreeMix infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (TreeMix infile)
#' @export TreeMix infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a TreeMix formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2treemix(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "TreeMix_infile.txt")
#' vcf2treemix(my_vcf, ind_pop, keepers, out_file = "TreeMix_infile.txt")
#' vcf2treemix(my_vcf, ind_pop, keepers)
#'

vcf2treemix <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE, 
                          out_file = "treemix_infile.txt")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf_list <- lapply(keep_pop, function(x) {
    vcf[, c(TRUE, x == ind_pop)]
  })
  names(vcf_list) <- keep_pop
  
  for (i in 1:length(vcf_list)) {
    gt <- extract.gt(vcf_list[[i]], return.alleles = F, convertNA = T) #convertNA not working here
    allele1 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 1, 1)
    })
    allele2 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 3, 3)
    })
    
    REF <- apply(cbind(allele1, allele2), 1, function(x){sum(x == 0, na.rm = TRUE)})
    ALT <- apply(cbind(allele1, allele2), 1, function(x){sum(x == 1, na.rm = TRUE)})
    GEN <- apply(cbind(REF, ALT), 1, function(x){glue::glue_collapse(x, sep = ",")})
    
    if (i == 1) {
      df <- GEN
    } else {
      df <- cbind(df, GEN)
    }
  }
  
  colnames(df) <- keep_pop
  utils::write.table(df, file = out_file, quote = FALSE, sep = " ", col.names = TRUE, row.names = FALSE)

  return(invisible(NULL))
}


################################
#' @title vcf2eigenstrat
#' @description converts vcfR format data to Eigenstrat infiles
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek December 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Eigenstrat infiles (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Eigenstrat infiles)
#' @export Eigenstrat infiles (3) of SNP positions, individuals, genotype matrix 
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Eigenstrat formatted input files
#' When list of sexes is not provided, lists all individuals as unknown
#' When relative position on chromosome (cM distance or similar) is not provides, list as 0
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2eigenstrat(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, sex = list_of_sex, rel_pos = marker_cM_map, inc_missing = TRUE, out_file = "Eigenstrat")
#' vcf2eigenstrat(my_vcf, ind_pop, keepers, out_file = "Eigenstrat")
#' vcf2eigenstrat(my_vcf, ind_pop, keepers)
#'

vcf2eigenstrat <- function (vcf, ind_pop, keep_pop, sex = "U", rel_pos = 0, inc_missing = TRUE, 
                         out_file = "eigenstrat_infile")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_sub_pops(vcf, ind_pop, keep_pop)
  gt <- extract.gt(vcf2, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "9"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "2"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "1"
  
  df <- apply(gt, 1, function(x){glue::glue_collapse(x)})
  utils::write.table(df, file = paste(out_file, ".geno", sep = ""), quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
  
  df <- cbind(colnames(gt), sex, as.character(ind_pop))
  utils::write.table(df, file = paste(out_file, ".ind", sep = ""), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  #rel_pos <- getPOS(vcf2)/1000
  df <- cbind(getID(vcf2), getCHROM(vcf2), rel_pos, getPOS(vcf2), getREF(vcf2), getALT(vcf2))
  utils::write.table(df, file = paste(out_file, ".snp", sep = ""), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  return(invisible(NULL))
}


################################
#' @title vcf2fineRadStructure
#' @description converts vcfR format data to fineRadStructure infile
#' @description in part based on vcfR2genepop function
#' @author Tomas Hrbek July 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in fineRadStructure infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (fineRadStructure infile)
#' @export fineRadStructure infile of alleles
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a fineRadStructure formatted input file
#' The function expects multiallelic loci
#'
#' @example
#' vcf2fineRadStructure(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "fineRadStructure_infile.txt")
#' vcf2fineRadStructure(my_vcf, ind_pop, keepers, out_file = "fineRadStructure_infile.txt")
#' vcf2fineRadStructure(my_vcf, ind_pop, keepers)
#'

vcf2fineRadStructure <- function (vcf, ind_pop, keep_pop, inc_missing = TRUE, 
                                  out_file = "fineRadStructure_infile.txt")
{
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums((is.na(gt))), ]
  }
  vcf2 <- vcf_sub_pops(vcf, ind_pop, keep_pop)
  # extract chromosome names
  chrom <- stringr::str_extract(getID(vcf2), "^[^_]*") %>%
    unique()
  # get list of SNPs by chromosomes
  snp_list <- lapply(chrom, function(x) {
    vcf2[x == stringr::str_extract(getID(vcf2), "^[^_]*"), ]
  })
  
  gt <- extract.gt(vcf2, return.alleles = F, convertNA = T)
  # write IDs to file
  write.table(gt[0,], quote = FALSE, sep = "\t", file = out_file)
  
  for (i in 1:length(snp_list)) {
    gt <- extract.gt(snp_list[[i]], return.alleles = T, convertNA = T) #convertNA not working here
    
    allele1 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 1, 1)
    })
    rownames(allele1) <- NULL
    allele1 <- t(allele1)
    allele1[allele1 == "."] <- NA
    allele1a <- apply(allele1, MARGIN = 1, function(x) {
      glue::glue_collapse(x)
    })
    allele2 <- apply(gt, MARGIN = 2, function(x) {
      substr(x, 3, 3)
    })
    rownames(allele2) <- NULL
    allele2 <- t(allele2)
    allele2[allele2 == "."] <- NA
    allele2a <- apply(allele2, MARGIN = 1, function(x) {
      glue::glue_collapse(x)
    })
    allele <- data.frame(allele1a, allele2a)
    allele$frs <- if_else(allele$allele1a == allele$allele2a, 
                          allele$allele1a, paste0(allele$allele1a, "/", allele$allele2a))
    allele$frs[is.na(allele$frs)] <- ""
    
    write.table(t(allele$frs), quote = FALSE, col.names = FALSE, row.names = FALSE, 
                sep = "\t", file = out_file, append = TRUE)
  }
  
  return(invisible(NULL))
}


################################
#' @title vcf2snapp
#' @description converts vcfR format data to SNAPP (Nexus) format infile
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in SNAPP infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (SNAPP infile)
#' @export SNAPP infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a SNAPP (Nexus) formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2snapp(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "SNAPP_infile.nex")
#' vcf2snapp(my_vcf, ind_pop, keepers, out_file = "SNAPP_infile.nex")
#' vcf2snapp(my_vcf, ind_pop, keepers)
#'

vcf2snapp <-function (vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "snapp.nex") 
{
  #    if (!require(ape)) {
  #        install.packages("ape")
  #    }
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_sub_pops(vcf, ind_pop, keep_pop)
  gt <- extract.gt(vcf2, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "?"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "1"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "2"
  gt <- t(gt)
  
  ape::write.nexus.data(gt, out_file, format = "standard", interleaved = FALSE)
  
  # fix symbols in nexus so that snapp.xml file is correctly formated
  fix_symbols <- paste(c("sed -i 's/symbols=\"0123456789\"/symbols=\"012\"/'", out_file), collapse = " ")
  system(fix_symbols)
  
  return(invisible(NULL))
}


################################
#' @title vcf2nexus
#' @description converts vcfR format data to Nexus format infile
#' @author Tomas Hrbek February 2021
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Nexus infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Nexus infile)
#' @export Nexus infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Nexus formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2nexus(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, inc_missing = TRUE, out_file = "Nexus_infile.nex")
#' vcf2nexus(my_vcf, ind_pop, keepers, out_file = "Nexus_infile.nex")
#' vcf2nexus(my_vcf, ind_pop, keepers)
#'

vcf2nexus <-function (vcf, ind_pop, keep_pop, inc_missing = TRUE, out_file = "nexus.nex") 
{
  #    if (!require(ape)) {
  #        install.packages("ape")
  #    }
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_sub_pops(vcf, ind_pop, keep_pop)
  gt <- extract.gt(vcf2, return.alleles = T, convertNA = T)
  gt[gt == "."] <- "?"
  gt[gt == "A/A" | gt == "A|A"] <- "A"
  gt[gt == "A/C" | gt == "C/A" | gt == "A|C" | gt == "C|A"] <- "M"
  gt[gt == "A/G" | gt == "G/A" | gt == "A|G" | gt == "G|A"] <- "R"
  gt[gt == "A/T" | gt == "T/A" | gt == "A|T" | gt == "T|A"] <- "W"
  gt[gt == "C/C" | gt == "C|C"] <- "C"
  gt[gt == "C/G" | gt == "G/C" | gt == "C|G" | gt == "G|C"] <- "S"
  gt[gt == "C/T" | gt == "T/C" | gt == "C|T" | gt == "T|C"] <- "Y"
  gt[gt == "G/G" | gt == "G|G"] <- "G"
  gt[gt == "G/T" | gt == "T/G" | gt == "G|T" | gt == "T|G"] <- "K"
  gt[gt == "T/T" | gt == "T|T"] <- "T"
  gt <- t(gt)
  
  ape::write.nexus.data(gt, out_file, format = "DNA", interleaved = FALSE)
  
  return(invisible(NULL))
}


################################
#' @title vcf2fasta
#' @description converts vcfR format data to Fasta format infile
#' @author Tomas Hrbek August 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Fasta infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Fasta infile)
#' @export Fasta infile of SNPs
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Fasta formatted input file
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @example
#' vcf2fasta(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, interleaved = FALSE, inc_missing = TRUE, out_file = "Fasta_infile.fas")
#' vcf2fasta(my_vcf, ind_pop, keepers, out_file = "Fasta_infile.fas")
#' vcf2fasta(my_vcf, ind_pop, keepers)
#'

vcf2fasta <-function (vcf, ind_pop, keep_pop, interleaved = FALSE, inc_missing = TRUE, out_file = "fasta.fas") 
{
  #    if (!require(ape)) {
  #        install.packages("ape")
  #    }
  if (class(vcf) != "vcfR") {
    stop(paste("Expecting an object of class vcfR, received a", 
               class(vcf), "instead"))
  }
  if (class(ind_pop) != "factor" | class(keep_pop) != "factor") {
    stop(paste("Expecting population vector, received a", 
               class(ind_pop), "and", class(keep_pop), "instead"))
  }
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  if (inc_missing == FALSE) {
    gt <- extract.gt(vcf, convertNA = T)
    vcf <- vcf[!rowSums(is.na(gt)), ]
  }
  vcf2 <- vcf_sub_pops(vcf, ind_pop, keep_pop)
  gt <- extract.gt(vcf2, return.alleles = T, convertNA = T)
  gt[gt == "."] <- "?"
  gt[gt == "A/A" | gt == "A|A"] <- "A"
  gt[gt == "A/C" | gt == "C/A" | gt == "A|C" | gt == "C|A"] <- "M"
  gt[gt == "A/G" | gt == "G/A" | gt == "A|G" | gt == "G|A"] <- "R"
  gt[gt == "A/T" | gt == "T/A" | gt == "A|T" | gt == "T|A"] <- "W"
  gt[gt == "C/C" | gt == "C|C"] <- "C"
  gt[gt == "C/G" | gt == "G/C" | gt == "C|G" | gt == "G|C"] <- "S"
  gt[gt == "C/T" | gt == "T/C" | gt == "C|T" | gt == "T|C"] <- "Y"
  gt[gt == "G/G" | gt == "G|G"] <- "G"
  gt[gt == "G/T" | gt == "T/G" | gt == "G|T" | gt == "T|G"] <- "K"
  gt[gt == "T/T" | gt == "T|T"] <- "T"
  gt <- t(gt)
  
  ifelse(interleaved == FALSE, nbcol <- -1, nbcol <- 6)
  
  ape::write.dna(gt, out_file, format = "fasta", nbcol = nbcol, colsep = "")
  
  return(invisible(NULL))
}


################################
#' @title vcf_sub_pops
#' @description subsets vcfR format data by population
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include/exclude in vcf infile (factor)
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by population (groups of individuals
#' assigned to one or more populations), returning new vcfR object
#'
#' @example
#' vcf_sub_pops(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, whitelist = TRUE)
#' vcf_sub_pops(vcf = my_vcf, ind_pop = ind_pop, keep_pop = non_keepers, whitelist = FALSE)
#' vcf_sub_pops(my_vcf, ind_pop, keepers)
#'

vcf_sub_pops <- function(vcf, ind_pop, keep_pop, whitelist = TRUE) {
  ids <- which(ind_pop %in% keep_pop)
  ids <- ids+1
  if (whitelist == TRUE){
    vcf <- vcf[, c(1,ids)] %>%
      vcf_filter_invariant()
  }
  else {
    vcf <- vcf[, -c(ids)] %>%
      vcf_filter_invariant()
  }
  
  return(vcf)
}


################################
#' @title vcf_sub_indivs
#' @description subsets vcfR format data by individuals
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param indiv -> individuals to retain/drop in vcf (factor)
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by individuals, returning new vcfR object
#'
#' @example
#' vcf_sub_indivs(vcf = my_vcf, indiv = indivs_to_keep, whitelist = TRUE)
#' vcf_sub_indivs(vcf = my_vcf, indiv = indivs_to_drop, whitelist = FALSE)
#' vcf_sub_indivs(my_vcf, indivs_to_keep)
#'

vcf_sub_indivs <- function(vcf, indiv, whitelist = TRUE) {
  # read all sample names in vcf
  vcf_names <- colnames(vcf@gt)[-1]
  
  # allow empty indiv list - keep all individuals
  if (length(indiv) == 0){
    indiv <- vcf_names
    whitelist <- TRUE
  }
  
  ids <- which(vcf_names %in% indiv)
  ids <- ids+1
  if (whitelist == TRUE){
    vcf <- vcf[, c(1,ids)] %>%
      vcf_filter_invariant()
  }
  else {
    vcf <- vcf[, -c(ids)] %>%
      vcf_filter_invariant()
  }
  
  return(vcf)
}


################################
#' @title vcf_filter_missingness
#' @description subsets vcfR format data by % missing data
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @param miss_p -> max missing data per locus as decimal (numeric)
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by % missing data, returning new vcfR object
#'
#' @example
#' vcf_filter_missingness(vcf = my_vcf, miss_p = miss_p)
#' vcf_filter_missingness(my_vcf, miss_p)
#'

vcf_filter_missingness <- function(vcf, miss_p) {
  gt <- extract.gt(vcf, convertNA = T)
  # get number of samples in vcf
  n_samples <- ncol(gt)
  
  # keep only those loci with < % missing data
  vcf <- vcf[rowSums(is.na(gt)) < floor(n_samples*miss_p),]
  
  # print VCF matrix completeness
  vcf1 <- vcf_filter_oneSNP(vcf)
  gt <- extract.gt(vcf1, convertNA = T)
  missing_p <- sum(is.na(gt)) / length(gt)
  print(paste("final % missing data in VCF is", round(missing_p*100, 2), "%", sep = " "))
  
  return(vcf)
}


################################
#' @title vcf_filter_missing_indivs
#' @description remove from vcfR format data indivs with >% missing data
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @param miss_p -> max missing data per locus as decimal (numeric)
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by % missing data within 
#' and individual, returning new vcfR object
#'
#' @example
#' vcf_filter_missing_indivs(vcf = my_vcf, miss_p = miss_p)
#' vcf_filter_missing_indivs(my_vcf, miss_p)
#'

vcf_filter_missing_indivs <- function(vcf, miss_p) {
  gt <- extract.gt(vcf, convertNA = T)
  # get number of snps in vcf
  n_snps <- nrow(gt)
  
  # report which samples above threshold
  samples_removed <- colnames(gt)[colSums(is.na(gt)) > floor(n_snps*miss_p)]
  cat(paste("removed samples are:", samples_removed, "\n", sep = " "))
  
  # keep only those loci with < % missing data
  vcf <- vcf[, c(TRUE, colSums(is.na(gt)) < floor(n_snps*miss_p))] %>%
    vcf_filter_invariant()
  
  # print VCF matrix completeness
  vcf1 <- vcf_filter_oneSNP(vcf)
  gt <- extract.gt(vcf1, convertNA = T)
  missing_p <- sum(is.na(gt)) / length(gt)
  print(paste("final % missing data in VCF is", round(missing_p*100, 2), "%", sep = " "))
  
  return(vcf)
}


################################
#' @title vcf_filter_multiSNP
#' @description subsets vcfR format data keeping only loci with 2+ SNPs
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object keeping only loci with between min and max # of SNPs
#' per locus, returning new vcfR object
#' default min = 2 and max = 5 SNPs per locus
#' Recommended as input for fineRADstructure analyses
#'
#' @example
#' vcf_filter_multiSNP(vcf = my_vcf, minS = 2, maxS = 5)
#' vcf_filter_multiSNP(vcf = my_vcf)
#' vcf_filter_multiSNP(my_vcf)
#'

vcf_filter_multiSNP <- function(vcf, minS = 2, maxS = 5) {
  # read all loci names in vcf
  chrom <- getCHROM(vcf)
  id <-  getID(vcf)
  chrom_id <- cbind(chrom, id) %>%
    as_tibble()
  
  keeper <- chrom_id %>%
    count(chrom) %>%
    filter(n >= minS & n <= maxS) %>%
    select(-n) %>%
    as.matrix() %>%
    as.character()
  
  # keep only those loci with between minS and maxS SNPs
  vcf <- vcf[chrom %in% keeper, ]
  
  # keep only those loci with 2+ SNPs
  # vcf <- vcf[chrom %in% unique(chrom[duplicated(chrom)]),]
  
  return(vcf)
}


################################
#' @title vcf_filter_oneSNP
#' @description subsets vcfR format data keeping only 1 SNP per locus
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object keeping only 1 SNP per locus, returning new vcfR object
#' The first SNP independent of quality is taken (may mofify this in future)
#'
#' @example
#' vcf_filter_oneSNP(vcf = my_vcf)
#' vcf_filter_oneSNP(my_vcf)
#'

vcf_filter_oneSNP <- function(vcf) {
  # read all loci names in vcf
  chrom <- getCHROM(vcf)
  
  # keep only those loci with 1 SNP
  vcf <- vcf[!duplicated(chrom),]
  
  return(vcf)
}


################################
#' @title vcf_filter_invariant
#' @description remove invariant loci from vcfR format data
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function removes invariant loci from the vcfR object
#' This might be desirable after subsetting the vcf by individuals
#'
#' @example
#' vcf_filter_invariant(vcf = my_vcf)
#' vcf_filter_invariant(my_vcf)
#'

vcf_filter_invariant <- function(vcf) {
  gt <- extract.gt(vcf, convertNA = T)
  # remove invariant loci
  y <- vector(length = nrow(gt))
  for(i in 1:length(y)) {
    # keep if num unique loci > 1
    y[i] <- length(unique(na.omit(gt[i,]))) > 1
  }
  vcf <- vcf[y,]
  
  return(vcf)
}


################################
#' @title vcf_filter_coverage
#' @description remove genotypes below coverage threshold from vcfR format data
#' @author Tomas Hrbek November 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function removes genotypes below coverage threshold from the vcfR object
#'
#' @example
#' vcf_filter_coverage(vcf = my_vcf, cover = 10)
#' vcf_filter_coverage(my_vcf)
#'

vcf_filter_coverage <- function(vcf, cover = 10) {
  dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  vcf@gt[,-1][dp < cover] <- "./.:0:.,.,.:0,0:0,0"
  
  return(vcf)
}


################################
#' @title vcf_filter_quality
#' @description remove loci below quality threshold from vcfR format data
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci below quality threshold from the vcfR object
#'
#' @example
#' vcf_filter_quality(vcf = my_vcf, qual = 20)
#' vcf_filter_quality(my_vcf)
#'

vcf_filter_quality <- function(vcf, qual = 20) {
  if(any(is.na(getQUAL(vcf)))){
    print("No quality information in VCF; keeping VCF as is")
  } else {
    # keep only those loci with minimum quality
    vcf <- vcf[getQUAL(vcf) >= qual,]
  }
  
  return(vcf)
}


################################
#' @title vcf_filter_rank
#' @description remove loci below rank threshold from vcfR format data
#' @description rank is calculated as sqrt(chi-sqr/n) of allele read counts
#' @description used for paralog detection - very low rank values (<0.4)
#' @description rank is calculated in DiscoSNP, registered as Pk in INFO
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci below rank threshold from the vcfR object
#'
#' @example
#' vcf_filter_quality(vcf = my_vcf, rank = .4)
#' vcf_filter_quality(my_vcf)
#'


vcf_filter_rank <- function(vcf, rank = .4) {
  snp_rank <- stringr::str_extract(vcf@fix[,8], "Rk=[0-9|.]*") %>%
    stringr::str_extract("[^Rk=]+")
  # keep only those loci with minimum rank
  vcf <- vcf[snp_rank >= rank,]
  
  return(vcf)
}


################################
#' @title vcf_filter_hets
#' @description remove loci above het threshold from vcfR format data
#' @description high hets are indicative of paralogs
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci above het threshold from the vcfR object
#'
#' @example
#' vcf_filter_hets(vcf = my_vcf, hets = .4)
#' vcf_filter_hets(my_vcf)
#'


vcf_filter_hets <- function(vcf, hets = .4) {
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  gt <- extract.gt(vcf, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "NA"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "1"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "2"
  
  homo <- apply(gt, MARGIN = 1, function(x) {
    sum(x == 0 | x == 1, na.rm = TRUE)
  })
  hetero <- apply(gt, MARGIN = 1, function(x) {
    sum(x == 2, na.rm = TRUE)
  })
  
  vcf <- vcf[hetero/(homo+hetero) <= hets,]
  
  return(vcf)
}


################################
#' @title vcf_filter_maf
#' @description remove loci below MAF threshold from vcfR format data
#' @author Tomas Hrbek October 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function removes loci below MAF threshold from the vcfR object
#'
#' @example
#' vcf_filter_maf(vcf = my_vcf, maf = .05)
#' vcf_filter_maf(my_vcf)
#'


vcf_filter_maf <- function(vcf, maf = .05) {
  # remove NAs from a vector
  rm_na <- function(x) {
    y <- x[!is.na(x)]
    return(y)
  }
  
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf), ]
  gt <- extract.gt(vcf, return.alleles = F, convertNA = T)
  gt[is.na(gt)] <- "NA"
  gt[gt == "0/0" | gt == "0|0"] <- "0"
  gt[gt == "1/1" | gt == "1|1"] <- "1"
  gt[gt == "0/1" | gt == "0|1" | gt == "1/0" | gt == "1|0"] <- "2"
  
  alleles_all <- apply(gt, MARGIN = 1, function(x) {
    2*length(rm_na(x))
  })
  allele_rare <- apply(gt, MARGIN = 1, function(x) {
    sum(x == 2, na.rm = TRUE) + 2*sum(x == 1, na.rm = TRUE)
  })
  
  vcf <- vcf[allele_rare/alleles_all >= maf,]
  
  return(vcf)
}
