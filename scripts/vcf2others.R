#!/usr/bin/env Rscript

# adapted vcfR2migrate function to allow for inclusion of missing data in Migrate output
# modifications on line 31
# Tomas Hrbek August 2020

vcf2migrate <- function (vcf, pop, in_pop, out_file = "MigrateN_infile.txt", 
    method = c("N", "H")) 
{
    method <- match.arg(method, c("N", "H"), several.ok = FALSE)
    if (class(vcf) != "vcfR") {
        stop(paste("Expecting an object of class vcfR, received a", 
            class(vcf), "instead"))
    }
    if (class(pop) != "factor") {
        stop(paste("Expecting population vector, received a", 
            class(pop), "instead"))
    }
    vcf <- extract.indels(vcf, return.indels = F)
    vcf <- vcf[is.biallelic(vcf), ]
    vcf_list <- lapply(in_pop, function(x) {
        vcf[, c(TRUE, x == pop)]
    })
    names(vcf_list) <- in_pop
    if (method == "N") {
        myHeader <- c("N", length(vcf_list), nrow(vcf_list[[1]]))
        pop_list <- vector(mode = "list", length = length(vcf_list))
        names(pop_list) <- names(vcf_list)
        for (i in 1:length(vcf_list)) {
            gt <- extract.gt(vcf_list[[i]], return.alleles = T)
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
