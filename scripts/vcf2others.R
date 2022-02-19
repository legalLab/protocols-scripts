#!/usr/bin/env Rscript

# functions to convert vcf data to other formats
# depends on vcfR (should be loaded since vcfR object [vcf file] is passed to these functions)
# depends on ape to generate nexus file (vcf2nexus, vcf2snapp)


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
    write(paste("#There are ", ncol(pop_list[[i]][[1]]), " SNPs", sep = ""), file = out_file, append = TRUE)
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
#' @param keep_pop -> population(s) of interest to include in Structure infile (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Structure infile)
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
        gt[gt == "A/A"] <- "0101"
        gt[gt == "A/C" | gt == "C/A"] <- "0102"
        gt[gt == "A/G" | gt == "G/A"] <- "0103"
        gt[gt == "A/T" | gt == "T/A"] <- "0104"
        gt[gt == "C/C"] <- "0202"
        gt[gt == "C/G" | gt == "G/C"] <- "0203"
        gt[gt == "C/T" | gt == "T/C"] <- "0204"
        gt[gt == "G/G"] <- "0303"
        gt[gt == "G/T" | gt == "T/G"] <- "0304"
        gt[gt == "T/T"] <- "0404"
        gt <- t(gt)
        write(paste("pop ", names(pop_list)[i], sep = ""), file = out_file, append = TRUE)
        write.table(cbind(sep = ',', gt), file = out_file, quote = FALSE, sep = " ", col.names = FALSE, append = TRUE)
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
    gt[gt == "0/0"] <- "0"
    gt[gt == "1/1"] <- "1"
    gt[gt == "0/1"] <- "2"
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
    gt[gt == "A/A"] <- "A"
    gt[gt == "A/C" | gt == "C/A"] <- "M"
    gt[gt == "A/G" | gt == "G/A"] <- "R"
    gt[gt == "A/T" | gt == "T/A"] <- "W"
    gt[gt == "C/C"] <- "C"
    gt[gt == "C/G" | gt == "G/C"] <- "S"
    gt[gt == "C/T" | gt == "T/C"] <- "Y"
    gt[gt == "G/G"] <- "G"
    gt[gt == "G/T" | gt == "T/G"] <- "K"
    gt[gt == "T/T"] <- "T"
    gt <- t(gt)
    
    ape::write.nexus.data(gt, out_file, format = "DNA", interleaved = FALSE)
    
    return(invisible(NULL))
}


################################
#' @title vcf_sub_pops
#' @description subsets vcfR format data by population
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in vcf infile (factor)
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by population (groups of individuals
#' assigned to one or more populations), returning new vcfR object
#'
#' @example
#' vcf_sub_pops(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers)
#' vcf_sub_pops(my_vcf, ind_pop, keepers)
#'

vcf_sub_pops <- function(vcf, ind_pop, keep_pop) {
    ids <- which(ind_pop %in% keep_pop)
    ids <- ids+1
    vcf <- vcf[, c(1,ids)]
    
    return(vcf)
}


################################
#' @title vcf_sub_indivs
#' @description subsets vcfR format data by individuals
#' @author Tomas Hrbek December 2020
#'
#' @param vcf -> vcfR object
#' @param indiv -> individuals to retain in vcf (factor)
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object by individuals, returning new vcfR object
#'
#' @example
#' vcf_sub_indivs(vcf = my_vcf, indiv = indivs_to_keep)
#' vcf_sub_indivs(my_vcf, indivs_to_keep)
#'

vcf_sub_indivs <- function(vcf, indiv) {
    # read all sample names in vcf
    vcf_names <- colnames(vcf@gt)[-1]
    
    ids <- which(vcf_names %in% indiv)
    ids <- ids+1
    vcf <- vcf[, c(1,ids)]
    
    return(vcf)
}


################################
#' @title filter_missingness
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
#' filter_missingness(vcf = my_vcf, miss_p = miss_p)
#' filter_missingness(my_vcf, miss_p)
#'

filter_missingness <- function(vcf, miss_p) {
    gt <- extract.gt(vcf, convertNA = T)
    # get number of samples in vcf
    n_samples <- ncol(vcf@gt) - 1
    
    # keep only those loci with < % missing data
    vcf <- vcf[rowSums(is.na(gt)) < floor(n_samples*miss_p), ]
    
    return(vcf)
}


################################
#' @title filter_multiSNP
#' @description subsets vcfR format data keeping only loci with 2+ SNPs
#' @author Tomas Hrbek February 2022
#'
#' @param vcf -> vcfR object
#' @export nothing
#' @return subsetted vcfR object
#'
#' @details
#' This function subsets the vcfR object keeping only loci with 2+ SNPs, returning new vcfR object
#' Recommended as input for fineRADstructure analyses
#'
#' @example
#' filter_multiSNP(vcf = my_vcf)
#' filter_multiSNP(my_vcf)
#'

filter_multiSNP <- function(vcf) {
    # read all loci names in vcf
    chrom <- getCHROM(vcf)
    
    # keep only those loci with 2+ SNPs
    vcf <- vcf[chrom %in% unique(chrom[duplicated(chrom)]) , ]
    
    return(vcf)
}
