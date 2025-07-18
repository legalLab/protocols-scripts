#!/usr/bin/env Rscript

# collection of helper functions
# Tomas Hrbek starting July 2019


# harmonic mean
harmonic_mean <- function(x) {
  x <- x[!is.na(x)]
  y <- 1/mean(1/x)
  return(y)
}

# geometric mean
geometric_mean <- function(x) {
  x <- x[!is.na(x)]
  y <- prod(x)^(1/length(x))
  return(y)
}

# modal density value
density_mode <- function(x) {
  x <- x[!is.na(x)]
  td <- density(x)
  max_dens <- which.max(td$y)
  y <- td$x[max_dens]
  return(y)
}

# mode of discrete series
# mode function from https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
discrete_mode <- function(x) {
  x <- x[!is.na(x)]
  ux <- unique(x)
  tab <- tabulate(match(x, ux)); 
  y <- ux[tab == max(tab)]
  return(y)
}

# remove NAs from a vector
rm_na <- function(x) {
  y <- x[!is.na(x)]
  return(y)
}

# function from https://github.com/tidyverse/tidyr/issues/279
explode <- function(df, w, .id = NULL) {
  w <- pull(mutate(df, `_weight` = !!enquo(w)))
  if (!is.numeric(w)) {
    stop("`w` must evalute to a numeric vector", call. = FALSE)
  }
  df <- df[rep(seq_len(nrow(df)), w), , drop = FALSE]
  if (!is.null(.id)) {
    df[[.id]] <- sequence(w)
  }
  return(df)
}

# generate table of samples and barcode sequences from table of samples and barcode names
recode_by_lookup <- function(df, lookup, ngs = "novogene", ...) {
  if(!is.data.frame(lookup) && !is.matrix(lookup)) {
    stop("The lookup table must be a data.frame or matrix")
  }
  if(ncol(lookup) != 2) {
    stop("The lookup table must have two colums")
  }
  if(length(unique(lookup[, 1])) < nrow(lookup) ) {
    stop("Index names in lookup table are not unique")
  }
  if(length(unique(lookup[, 2])) < nrow(lookup) ) {
    stop("Index sequences in lookup table are not unique")
  }
  if(!is.data.frame(df)) {
    stop("The dataframe must be a tibble or data.frame")
  }
  if(ncol(df) != 3 && ngs == "azenta") {
    stop("The Illumina dataframe must have three colums; q5 q7 id")
  }
  else if (ncol(df) != 3 && ngs == "novogene") {
    stop("The Illumina dataframe must have three colums; q5 q7 id")
  }
  else if (ncol(df) != 2 && ngs == "ion") {
    stop("The IonTorrent dataframe must have two colums; bc id")
  }
  else if (ncol(df) != 3 && ngs == "ion_seq") {
    stop("The IonTorrent dataframe must have three colums; id A P1")
  }
  else if (ncol(df) != 4 && ngs == "ref_seq") {
    stop("The Illumina dataframe must have four colums; q5 q7 pos id")
  }
  else if (ncol(df) != 4 && ngs == "edna") {
    stop("The Illumina dataframe must have four colums; q5 q7 pos id")
  }
  else if (ncol(df) != 4 && ngs == "micro_seq") {
    stop("The Illumina dataframe must have four colums; q5 q7 pos id")
  }
  
  # convert tibble to dataframe (otherwise df operations do not work)
  df <- as.data.frame(df)
  # remove missing samples
  df <- df[!is.na(df$id), ]
  
  if(ngs == "azenta") {
    # lookup q5
    df$Q5n <- df$Q5
    df$Q5 <- lookup[match(df$Q5, lookup[, 1]), 2]
    # lookup q7
    df$Q7n <- df$Q7
    df$Q7 <- lookup[match(df$Q7, lookup[, 1]), 2]
    # reorder columns
    col_order <- c("id", "Q7n", "Q7", "Q5n", "Q5")
    df <- df[, col_order]
  }
  else if(ngs == "novogene") {
    # lookup q5
    df$Q5 <- lookup[match(df$Q5, lookup[, 1]), 2]
    # lookup q7
    df$Q7 <- lookup[match(df$Q7, lookup[, 1]), 2]
    # reorder columns
    col_order <- c("id", "Q7", "Q5")
    df <- df[, col_order]
  }
  else if(ngs == "ion") {
    # lookup ham
    df$bc <- lookup[match(df$bc, lookup[, 1]), 2]
    # reorder columns
    col_order <- c("id", "bc")
    df <- df[, col_order]
  }
  else if(ngs == "ion_seq") {
    # lookup A
    df$A <- lookup[match(df$A, lookup[, 1]), 2]
    # lookup P1
    df$P1 <- lookup[match(df$P1, lookup[, 1]), 2]
    # reorder columns
    col_order <- c("id", "A", "P1")
    df <- df[, col_order]
  }
  else if(ngs == "ref_seq") {
    # lookup q5
    df$Q5n <- df$Q5
    df$Q5 <- lookup[match(df$Q5, lookup[, 1]), 2]
    # lookup q7
    df$Q7n <- df$Q7
    df$Q7 <- lookup[match(df$Q7, lookup[, 1]), 2]
    # lookup Ham
    df$Hamn <- ham
    df$Ham <- lookup[match(df$Hamn, lookup[, 1]), 2]
    df$primerF <- primerF
    df$primerR <- primerR
    # reorder columns
    col_order <- c("Q7n", "Q7", "Q5n", "Q5", "Hamn", "Ham", "primerF", "primerR", "pos", "id")
    df <- df[, col_order]
  }
  else if(ngs == "edna") {
    # lookup q5
    df$Q5n <- df$Q5
    df$Q5 <- lookup[match(df$Q5, lookup[, 1]), 2]
    # lookup q7
    df$Q7n <- df$Q7
    df$Q7 <- lookup[match(df$Q7, lookup[, 1]), 2]
    df$primerF <- primerF
    df$primerR <- primerR
    # reorder columns
    col_order <- c("pos", "id", "primerF", "primerR")
    df <- df[, col_order]
  }
  else if(ngs == "micro_seq") {
    # lookup q5
    df$Q5n <- df$Q5
    df$Q5 <- lookup[match(df$Q5, lookup[, 1]), 2]
    # lookup q7
    df$Q7n <- df$Q7
    df$Q7 <- lookup[match(df$Q7, lookup[, 1]), 2]
    # subset primer list by taxon
    primers <- primers[primers[,1] == taxon,]
    # lookup micros, loop through primers
    for (i in 1:nrow(primers)) {
      df1 <- df
      df1$taxon <- primers[i,1]
      df1$primer <- primers[i,2]
      df1$primerF <- primers[i,3]
      df1$primerR <- primers[i,4]
      if (i == 1) {
        df2 <- df1
      } else {
        df2 <- add_row(df2, df1)
      }
    }
    # reorder columns
    col_order <- c("Q7n", "Q7", "Q5n", "Q5", "taxon", "primer", "primerF", "primerR", "pos", "id")
    df <- df2[, col_order]
  }
  return(df)
}

# process NCBI local blast+ output
analyze_edna <- function(infile, path, lookup_table, pid = .97, len = 99, n_reads = 4) {
  #check if file is not empty first - empty gz files are 61-62 bytes
  #try(if(file.size(paste0(path, infile)) <= 100L) stop(paste0("empty file: ", infile)))
  
  # check if file is not empty first - empty gz files are 61-62 bytes
  # make a fake file with 1 sample "bob"
  if(file.size(paste0(path, infile)) <= 100L) {
    blast_table <- tibble(id = "bob", n = 1, freq = 0, taxonomy = "Bob;Bob;Bob;Bob;Bob;Bob;Bob bob")
    return(blast_table)
    break
  }
  
  blast <- read.table(gzfile(paste0(path, infile)), header = FALSE, sep = "\t")
  colnames(blast) <- c('qacc', 'sacc', 'evalue', 'pident', 'mismatch', 'length', 'qstart', 'qend', 'sstart', 'send')
  
  blast <- as_tibble(blast) %>%
    distinct(qacc, .keep_all = TRUE) %>%
    mutate(pident = pident/100) %>%
    filter(pident >= pid & length > len)
  
  #blast$sacc <- str_extract(blast$sacc, regex("[a-zA-Z0-9_]+$"))
  blast_table <- blast %>%
    mutate(n = as.numeric(str_extract(.$qacc, pattern = "[0-9]+$"))) %>%
    group_by(sacc) %>%
    summarise(n = sum(n)) %>%
    mutate(freq = n/sum(n)) %>%
    filter(n > n_reads) %>%
    arrange(desc(n)) %>%
    rename(id = sacc) %>%
    left_join(lookup_table, by = 'id') %>%
    add_row(id = "bob", n = 1, freq = 0, taxonomy = "Bob;Bob;Bob;Bob;Bob;Bob;Bob bob")
  
  return(blast_table)
}
