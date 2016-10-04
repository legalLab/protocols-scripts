# R script to convert PyRAD .loci files into input for G-PhoCS

# load libs - need to make sure it's pointing to the correct package
require("ape", lib.loc="/home/rupert/Software/Rpackages")

# load up the sample names
tab <- sort(scan(file="samples_names.txt", what="character"))

# get names of all the locus fasta files
fl <- list.files(path=".", pattern="locus")

# read fasta
dfl <- lapply(fl, function(x) read.dna(file=x, format="fasta", comment.char="#"))

# prep for loop
np <- list()
dummies <- list()
all <- list()
sall <- list()
write(length(fl), file="gphocs.loci", append=FALSE)
# run loop to add the missing data and create new locus file
for (i in 1:length(fl)){#
    if (setequal(tab, labels(dfl[[i]])) == TRUE){#
        write(paste0(fl[i], " #"), file="gphocs.loci", append=TRUE)
        write.dna(dfl[[i]][sort(labels(dfl[[i]])), ], file="gphocs.loci", format="interleaved", append=TRUE, colw=9999)
        }#
    else {#
        np[[i]] <- setdiff(tab, labels(dfl[[i]]))
        dummies[[i]] <- matrix(data="n", nrow=length(np[[i]]), ncol=dim(dfl[[i]])[2])
        row.names(dummies[[i]]) <- np[[i]]
        dummies[[i]] <- as.DNAbin(dummies[[i]])
        all[[i]] <- rbind(dfl[[i]], dummies[[i]])
        sall[[i]] <- all[[i]][sort(labels(all[[i]])), ]
        write(paste0("\n", fl[i], " #"), file="gphocs.loci", append=TRUE)
        write.dna(sall[[i]], file="gphocs.loci", format="interleaved", append=TRUE, colw=9999)
        }#
}#
#
