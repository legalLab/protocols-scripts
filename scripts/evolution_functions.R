#!/usr/bin/env Rscript

# functions for tracking selection in a three allele Hemoglobin system
# based on Templeton, A.R., 2006. Population Genetics and Microevolutionary Theory. John Wiley & Sons, New York, NY.

################################
#' @title genotype_f
#' @description calculates genotypic frequencies from allelic frequencies
#' @author Tomas Hrbek January 2021
#'
#' @param f_all -> frequencies of alleles (collection of numeric)
#' @export nothing
#' @return f_gen -> frequencies of genotypes (collection of numeric)
#'
#' @details
#' This function takes frequencies of three alleles and calculates the expected
#' HWE genotypic frequencies (six genotype)
#' The alleles are A, S, C
#' The genotypes are AA, AS, AC, SS, SC, CC
#'
#' @example
#' my_genotypes <- genotype_f(c(.5, .2, .3))
#' my_alelles <- c(.5, 0, .5)
#' my_genotypes <- genotype_f(my_alelles)
#'

genotype_f <- function(f_all)
{
  # if any allele has frequency less then 0.001 (0.1%) make 0
  if (f_all[1] <= 0.001) f_all[1] <- 0
  if (f_all[2] <= 0.001) f_all[2] <- 0
  if (f_all[3] <= 0.001) f_all[3] <- 0
  
  # if allelic frequencies do not sum up to 1, adjust highest frequency allele
  if (sum(f_all) != 1) {
    f_all[f_all == max(f_all)] <- max(f_all) + 1-sum(f_all)
  }

  fA <- f_all[1]
  fS <- f_all[2]
  fC <- f_all[3]
  
  # genotypic frequencies
  f_gen <- c(1:6)
  
  f_gen[1] <- fA^2
  f_gen[2] <- 2*fA*fS
  f_gen[3] <- 2*fA*fC
  f_gen[4] <- fS^2
  f_gen[5] <- 2*fS*fC
  f_gen[6] <- fC^2
  
  return(f_gen)
}


################################
#' @title allele_f
#' @description calculates allelic frequencies from genotypic frequencies
#' @author Tomas Hrbek January 2021
#'
#' @param f_gen -> frequencies of genotypes (collection of numeric)
#' @param w_gen -> adaptive values of genotypes (collection of numeric)
#' @export nothing
#' @return f_all -> frequencies of alleles (collection of numeric)
#'
#' @details
#' This function takes frequencies of six genotypes and calculates the expected
#' HWE allelic frequencies (three alleles)
#' The alleles are A, S, C
#' The genotypes are AA, AS, AC, SS, SC, CC
#' The fitnesses are wAA, wAS, wAC, wSS, wSC, wCC
#'
#' @example
#' my_alleles <- allele_f(c(.1, .1, .1, .1, .1, .5))
#' my_genotypes <- c(.1, 0, .4, 0, 0, .5)
#' my_alleles <- allele_f(my_genotypes)
#'

allele_f <- function(f_gen, w_gen) {

  fAA <- f_gen[1]
  fAS <- f_gen[2]
  fAC <- f_gen[3]
  fSS <- f_gen[4]
  fSC <- f_gen[5]
  fCC <- f_gen[6]
  wAA <- w_gen[1]
  wAS <- w_gen[2]
  wAC <- w_gen[3]
  wSS <- w_gen[4]
  wSC <- w_gen[5]
  wCC <- w_gen[6]
  
  # adaptive value of an average phenotype
  w <- sum(f_gen*w_gen)
  
  # allelic frequencies
  f_all <- c(1:3)
  
  f_all[1] <- fAA*wAA/w + .5*fAS*wAS/w + .5*fAC*wAC/w
  f_all[2] <- fSS*wSS/w + .5*fAS*wAS/w + .5*fSC*wSC/w
  f_all[3] <- fCC*wCC/w + .5*fAC*wAC/w + .5*fSC*wSC/w
  
  # make sure once one allele is 1, all others are 0 (floating point precission issues)
  if(f_all[1] == 1) f_all[c(2,3)] <- 0
  if(f_all[2] == 1) f_all[c(1,3)] <- 0
  if(f_all[3] == 1) f_all[c(1,2)] <- 0
  
  return(f_all)
}


################################
#' @title average_excess
#' @description calculates the average excess of individual alleles
#' @author Tomas Hrbek January 2021
#'
#' @param f_all -> frequencies of alleles (collection of numeric)
#' @param w_gen -> adaptive values of genotypes (collection of numeric)
#' @export nothing
#' @return w -> population average fitness (numeric)
#' @return a_all -> average excess of alleles (collection of numeric)
#'
#' @details
#' This function takes frequencies of three alleles and the adaptive values of 
#' genotypes formed by these alleles to calculate the average excess of each allele
#' The alleles are A, S, C
#' The genotypes are AA, AS, AC, SS, SC, CC
#' The fitnesses are wAA, wAS, wAC, wSS, wSC, wCC
#'
#' @example
#' my_alleles <- allele_f(c(.1, .1, .1, .1, .1, .5))
#' my_genotypes <- c(.1, 0, .4, 0, 0, .5)
#' my_alleles <- allele_f(my_genotypes)
#'

average_excess <- function(f_gen, w_gen) {

  fAA <- f_gen[1]
  fAS <- f_gen[2]
  fAC <- f_gen[3]
  fSS <- f_gen[4]
  fSC <- f_gen[5]
  fCC <- f_gen[6]
  wAA <- w_gen[1]
  wAS <- w_gen[2]
  wAC <- w_gen[3]
  wSS <- w_gen[4]
  wSC <- w_gen[5]
  wCC <- w_gen[6]
  
  # adaptive value of an average phenotype
  w <- sum(f_gen*w_gen)
  
  # allelic frequencies
  fA <- fAA + .5*fAS + .5*fAC
  fS <- fSS + .5*fAS + .5*fSC
  fC <- fCC + .5*fAC + .5*fSC
  
  # average excess of alleles
  a_all <- c(1:3)
  
  if(fA != 0) a_all[1] <- fA*(wAA-w) + fS*(wAS-w) + fC*(wAC-w) else a_all[1] <- NA
  if(fS != 0) a_all[2] <- fS*(wSS-w) + fA*(wAS-w) + fC*(wSC-w) else a_all[2] <- NA
  if(fC != 0) a_all[3] <- fC*(wCC-w) + fA*(wAC-w) + fS*(wSC-w) else a_all[3] <- NA
  
  return(c(w, a_all))
}


################################
#' @title drift
#' @description simulates genetic drift over one generation
#' @author Tomas Hrbek January 2021
#'
#' @param f_all -> frequencies of alleles (collection of numeric)
#' @param Ne -> effective population size (integer)
#' @export nothing
#' @return f_all -> frequencies of alleles (collection of numeric)
#'
#' @details
#' This function generates of collection of Ne alleles given f_all frequencies
#' of individual alleles, and samples this collection with replacement Ne times
#' From the resampled collection, new allelic frequencies are calculated
#' The alleles are A, S, C
#'
#' @example
#' my_alelles <- drift(c(.5, .2, .3), 100)
#' N <- 100
#' my_alelles <- c(.5, 0, .5)
#' my_alelles <- drift(my_alelles, N)
#'

drift <- function(f_all, Ne) {
  
  fA <- f_all[1]
  fS <- f_all[2]
  fC <- f_all[3]
  
  alleles <- c(rep("A",round(fA*Ne)), rep("S",round(fS*Ne)), rep("C",round(fC*Ne)))
  alleles <- sample(alleles, Ne, replace=TRUE)
  
  f_all[1] <- length(alleles[alleles == "A"])/Ne
  f_all[2] <- length(alleles[alleles == "S"])/Ne
  f_all[3] <- length(alleles[alleles == "C"])/Ne
  
  return(f_all)
}


################################
#' @title evolve
#' @description simulates selection over one generation with or without drift
#' @author Tomas Hrbek January 2021
#'
#' @param genotypes -> frequencies of genotypes (collection of numeric)
#' @param fitness -> fitness of genotypes (collection of numeric)
#' @param gens -> number of generations (integer)
#' @param Ne -> effective population size (integer); default = 0
#' @param drift -> implement genetic drift (logical)
#' @export nothing
#' @return df -> dataframe with genotypic and allelic frequencies, genotype and population average fitnesses
#' 
#' @details
#' This function calculates genotypic and allelic frequencies across gens 
#' generations give selection and genetic drift
#' Calls genotype_f, allele_f and drift functions
#' The alleles are A, S, C
#' The genotypes are AA, AS, AC, SS, SC, CC
#' The fitnesses are wAA, wAS, wAC, wSS, wSC, wCC
#'
#' @example
#' genotypes <- c(.1, 0, .4, 0, 0, .5)
#' fitness <- c(.9, 1, .9, .2, .7, 1.3)
#' my_df <- select(genotypes, fitness, 100)
#' my_df <- select(genotypes, fitness, 100, 10000, drift = TRUE)
#'

evolve <- function(genotypes, fitness, gens, Ne = 0, drift = FALSE) {
#  if(sum(genotypes) != 1) {
#    stop(paste("Genotypic frequencies must sum up to 1, instead they sum to", 
#               sum(genotypes)))
#  }
  
  if(drift == TRUE & Ne == 0) {
    stop(paste("Enter Ne larger than the default '0'"))
  }
  fitness_equal <- c(1, 1, 1, 1, 1, 1)
  alleles <- allele_f(genotypes, fitness_equal)
  excess <- average_excess(genotypes, fitness)
  df <- c(round(genotypes, 3), round(fitness, 3), round(excess, 3), round(alleles, 3))
  names(df) <- c('fAA', 'fAS', 'fAC', 'fSS', 'fSC', 'fCC', 'wAA', 'wAS', 'wAC', 
                 'wSS', 'wSC', 'wCC', 'w', 'aA', 'aS', 'aC', 'fA', 'fS', 'fC')
  
  for(rep in 2:gens) {
    if(drift == TRUE) {
      alleles <- drift(alleles, Ne)
    }
    genotypes <- genotype_f(alleles)
    alleles <- allele_f(genotypes, fitness)
    excess <- average_excess(genotypes, fitness)
    df1 <- c(round(genotypes, 3), round(fitness, 3), round(excess, 3), round(alleles, 3))
    df <- rbind(df, df1)
  }
  rownames(df) <- c(1:gens)
  
  #return(df) # returns "matrix" "array"
  return(as.data.frame(df))
}



################################
# small helper function to generate simple allele frequency plots
plot_alleles <- function(x) {
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(xpd=T, mar=par()$mar+c(0,0,0,3))
  matplot(x$fA, type = 'l', ylim = c(0,1), xlab = 'generations', ylab = 'allelic frequencies')
  matlines(x$fS, type = 'l', col = 'red')
  matlines(x$fC, type = 'l', col = 'blue')
  legend(1.07*nrow(x), 1, c('fA', 'fS', 'fC'), cex = 0.8, col = c('black', 'red', 'blue'), lty=c(1, 1, 1))
  
  return(invisible(NULL)) #just to keep things honest
}

################################
# small helper function to generate simple genotype frequency plots
plot_genotypes <- function(x) {
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(xpd=T, mar=par()$mar+c(0,0,0,3))
  matplot(x$fAA, type = 'l', ylim = c(0,1), xlab = 'generations', ylab = 'genotypic frequencies')
  matlines(x$fAS, type = 'l', col = 'red')
  matlines(x$fAC, type = 'l', col = 'blue')
  matlines(x$fSS, type = 'l', col = 'green')
  matlines(x$fSC, type = 'l', col = 'yellow')
  matlines(x$fCC, type = 'l', col = 'orange')
  legend(1.07*nrow(x), 1, c('fAA', 'fAS', 'fAC', 'fSS', 'fSC', 'fCC'), cex = 0.8, 
         col = c('black', 'red', 'blue', 'green', 'yellow', 'orange'), lty=c(1, 1, 1, 1, 1, 1))
  
  return(invisible(NULL)) #just to keep things honest
}

################################
# small helper function to generate simple allele frequency plots
plot_fitness <- function(x) {
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(xpd=T, mar=par()$mar+c(0,0,0,3))
  matplot(x$w, type = 'l', ylim = c(0,max(x$w)*1.1), xlab = 'generations', ylab = 'fitness', col = 'green')
  legend(1.07*nrow(x), 1, 'w', cex = 0.8, col = 'green', lty = 1)
  
  return(invisible(NULL)) #just to keep things honest
}


################################
# small helper function to summarize outcomes (mean and sd, counts of outcomes)
summarize_evolve <- function(genotypes, fitness, gens, Ne = 0, drift = FALSE, reps = 1000) {
  tmp_matrix <- matrix(nrow = reps, ncol = 10)
  sum_matrix <- matrix(nrow = 2, ncol = 17)
  for(rep in 1:reps) {
    evol <- evolve(genotypes, fitness, gens, Ne, drift)
    tmp_matrix[rep,] <- as.matrix(evol[gens, c(1:6, 13, 17:19)])
  }
  sum_matrix[1,1:10] <- round(apply(tmp_matrix, 2, mean), 3)
  sum_matrix[2,1:10] <- round(apply(tmp_matrix, 2, sd), 3)
  sum_matrix[1,11:13] <- apply(tmp_matrix[,8:10], 2, function(x) sum(x == 1))
  sum_matrix[1,17] <- sum(apply(tmp_matrix[,c(8,9,10)], 1, function(x) sum(x != 0) == 3))
  sum_matrix[1,14] <- sum(apply(tmp_matrix[,c(8,9)], 1, function(x) sum(x != 0) == 2)) - sum_matrix[1,17]
  sum_matrix[1,15] <- sum(apply(tmp_matrix[,c(8,10)], 1, function(x) sum(x != 0) == 2)) - sum_matrix[1,17]
  sum_matrix[1,16] <- sum(apply(tmp_matrix[,c(9,10)], 1, function(x) sum(x != 0) == 2)) - sum_matrix[1,17]
  colnames(sum_matrix) <- c('fAA', 'fAS', 'fAC', 'fSS', 'fSC', 'fCC', 'w', 'fA', 'fS', 'fC', 'fixed_A', 'fixed_S', 'fixed_C', 'polymorphic_AS', 'polymorphic_AC', 'polymorphic_CS', 'polymorphic_ACS')
  rownames(sum_matrix) <- c('mean', 'sd')
  
  return(as.data.frame(sum_matrix))
}

