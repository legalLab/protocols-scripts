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
recode_by_lookup <- function (df, lookup, ngs = "illumina") {
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
  else if (ncol(df) != 3 && ngs == "illumina") {
    stop("The Illumina dataframe must have three colums; q5 q7 id")
  }
  else if (ncol(df) != 2 && ngs == "ion") {
    stop("The IonTorrent dataframe must have two colums; bc id")
  }
  else if (ncol(df) != 3 && ngs == "ion_seq") {
    stop("The IonTorrent dataframe must have three colums; id A P1")
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
  else if(ngs == "illumina") {
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
  
  return(df)
}

