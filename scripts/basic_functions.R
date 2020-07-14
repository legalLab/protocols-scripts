#!/usr/bin/env Rscript

# functions to calculate harmoninc and geometric means
# Tomas Hrbek July 2019


# function - harmonic mean
harmonicMean <- function(logLikelihoods) {
  1/mean(1/logLikelihoods)
}
# function geometric mean
geometricMean <- function(logLikelihoods) {
  prod(logLikelihoods)^(1/length(logLikelihoods))
}

# examples
