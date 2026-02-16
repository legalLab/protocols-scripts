#library(ggplot2)
#library(dplyr)
#library(tidyr)

simulate_drift <- function(n, p, n_gen, n_sim) 
{
  #' simulate_drift(n, p, n_gen, n_sim)  
  #' 
  #' Simulate genetic drift using a Wright-Fisher model.
  #' 
  #' n: number of diploid individuals in the population
  #' p: starting allele frequency
  #' n_gen: number of generations
  #' n_sim: number of times to run the simulation
  #' 
  #' Returns a dataframe with n_gen rows and n_sim columns. Each 
  #' row represents one geneneration and each column are the results
  #' from one simulation. The value in each cell is the allele 
  #' frequency in that row's generation, in that column's run of 
  #' the simulation.
  
  q <- 1-p
  n_chrom <- 2*n # number of chromosomes
  
  drift <- array(0, dim=c(n_gen, n_sim))
  drift[1,] <- n_chrom*p # initialize number of A1 alleles in first generation
  for (j in 1:n_sim) {
    for(i in 2:n_gen) {
      drift[i,j] <- rbinom(1, n_chrom, prob=drift[i-1,j]/n_chrom)
    }  
  }
  drift <- data.frame(drift/n_chrom)
  return(drift)
}

plot_drift <- function(drift, n)
{
  #' plot_drift(drift)
  #' 
  #' Plot the result of simulate_drift()
  #' 
  #' drift: the dataframe that is returned by simulate_drift()
  #' 
  #' Returns a plot of the dataframe with each simulation a different 
  #' color.
  
  n_gen <- nrow(drift)
  n_sim <- ncol(drift)
  
  df_long <- drift %>%
    mutate(gens = 1:nrow(.)) %>%
    pivot_longer(cols = -gens, names_to = "reps", values_to = "freq")
  plt <- ggplot(df_long, aes(x = gens, y = freq, group = reps, colour = reps)) + 
    geom_line() + 
    labs(title = paste0("Simulations of Genetic Drift, n = ", n), x = "Generation", y = "Allele Frequency", colour = "Simulations")
  return(plt)
}

gens_to_fix <- function(n, p, n_rep)
{
  #' generations_to_fix(n, p, n_rep)
  #' 
  #' Simulate drift using a Wright-Fisher model and record how many
  #' generations until an allele is fixed or lost and whether the 
  #' allele was fixed or lost.
  #' 
  #' n: number of diploid individuals in the population
  #' p: the starting allele frequency
  #' n_rep: the number of times to repeat the simulation
  #' 
  #' Returns a dataframe of length n_rep with two columns, 
  #' "Generations" and "Fixed". Each row is the result of one 
  #' simulation.  The "Generations" column is an integer that records
  #' the number of generations until the allele was fixed or lost. 
  #' The "Fixed" column is a logical, either TRUE or FALSE. It is 
  #' TRUE if the allele was fixed and FALSE if the allele was lost.
  
  n_chrom <- 2*n # number of chromosomes/alleles
  gen <- vector("integer", length = n_rep)
  fixed <- vector("logical", length = n_rep)

  for (i in c(1:n_rep)) {
    j <- 0
    a <- n_chrom*p
    while (a > 0 & a < n_chrom) {
      a <- rbinom(1, n_chrom, prob = a/n_chrom)
      j <- j+1
    }
    gen[i] <- j
    fixed[i] <- a
  }
  return(data.frame(gen = gen, fixed = fixed))
}

rate_of_loss_fix <- function(df, p)
{
  #' rate_of_loss(df)
  #' 
  #' Computes the rate of allelic loss 
  #' 
  #' df: the dataframe returned by generations_to_fix()
  #' 
  #' Returns the rate of allelic loss (alleles per generation)
  
  p1 <- ifelse(p < .5, p, 1-p)
  r <- (1/mean(df$gen))/p1
  return(r)
}

plot_histogram <- function(df, n){
  # transform vector to tibble and add information
  df <- df %>%
    as_tibble()
  
  # plot a histogram
  plt <- ggplot(df, aes(x = gen)) +
    geom_histogram() +
    ggtitle(paste0("Histogram of Fixation Times, n = ", n)) +
    xlab("Time to Fixation")
  
  return(plt)
}

plot_log10_histogram <- function(df, n){
  # transform vector to tibble and add information
  df <- df %>%
    as_tibble()
  
  # plot a histogram
  plt <- ggplot(df, aes(x = gen)) +
    geom_histogram() +
    scale_x_log10() +
    ggtitle(paste0("Histogram of Fixation Times, n = ", n)) +
    xlab("Log10 time to Fixation")
  
  return(plt)
}

plot_boxplot <- function(df, n){
  # transform vector to tibble and add information
  df <- df %>%
    as_tibble() %>%
    mutate(run = 'run')
  
  # plot a boxplot
  plt <- ggplot(df, aes(x = gen, y = run)) +
    geom_boxplot() + 
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    ggtitle(paste0("Histogram of Fixation Times, n = ", n)) +
    xlab("Time to Fixation")
  
  return(plt)
}

plot_log10_boxplot <- function(df, n){
  # transform vector to tibble and add information
  df <- df %>%
    as_tibble() %>%
    mutate(run = 'run')
  
  # plot a boxplot
  plt <- ggplot(df, aes(x = gen, y = run)) +
    geom_boxplot() + 
    scale_x_log10() +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    ggtitle(paste0("Histogram of Fixation Times, n = ", n)) +
    xlab("Log10 time to Fixation")
  
  return(plt)
}
