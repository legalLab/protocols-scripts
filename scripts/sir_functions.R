#!/usr/bin/env Rscript

#------------------------------
# deterministic SIR function
sir_det <- function(N, t, b, a){
  sir_m <- matrix(nrow = t, ncol = 3) #make empty matrix
  
  colnames(sir_m) <- c("S", "I", "R") # assign names to columns
  rownames(sir_m) <- c(1:t)
  
  #initialize matrix with values of first iteration
  sir_m[1,1] <- N-1 # S = susceptible
  sir_m[1,2] <- 1 # I = infected
  sir_m[1,3] <- 0 # R = recovered
  
  #b is rate of infection
  #a is rate of recovery
  
  for (i in 2:t){
    if (sir_m[i-1,1] <= b*sir_m[i-1,1]*sir_m[i-1,2]){ #if # newly infected > susceptible (WHY *sir_m[i-1,3])
      # S = 0 if # newly infected individuals > S in the previous iteration
      sir_m[i,1] <- 0
      # I = I in the previous timestep + S in the previous iteration - newly recovered individuals in the previous iteration
      sir_m[i,2] <- sir_m[i-1,2] + sir_m[i-1,1] - a*sir_m[i-1,2]
    }
    else{
      # S = S in the previous iteration - newly infected individuals in the previous iteration
      sir_m[i,1] <- sir_m[i-1,1] - b*sir_m[i-1,1]*sir_m[i-1,2]
      # I = I in the previous timestep + newly infected individuals in the previous iteration - newly recovered individuals in the previous iteration
      sir_m[i,2] <- sir_m[i-1,2] + b*sir_m[i-1,1]*sir_m[i-1,2] - a*sir_m[i-1,2]
    }
    # R = R in the previous iteration + newly recovered individuals in the previous iteration
    sir_m[i,3] <- sir_m[i-1,3] + a*sir_m[i-1,2]
  }
  return(as.data.frame(sir_m))
}


#------------------------------
# SIR function that does the actual stochastic simulations
SIR <- function(n, t, b, a) {
  # n - number of individuals
  # t - length (iteration/time) of run
  # b - infection rate
  # a - recovery rate
  
  # create and initialize matrix
  dm <- matrix(nrow = t, ncol = n)
  dm[1,1] <- "I"  # sets the first infected individuals
  dm[1,2:n] <- "S" # all others are susceptible
  
  # loop through time t (rows)
  for (i in 2:t) {
    # count number of infected individuals
    inf <- sum(dm[i-1,] == "I")
    # loop through individuals (columns)
    for (j in 1:n) {
      # if individual is susceptible can become infected or stays susceptible
      if (dm[i-1,j] == "S") {
        if (runif(1) <= (1-(1-b)^inf)) {
          dm[i,j] <- "I"
        }
        else {dm[i,j] <- dm[i-1,j] # S
        }
      }
      # else if individual is infected can become removed or stays infected
      else if (dm[i-1,j] == "I") {
        if(runif(1) <= a) {
          dm[i,j] <- "R"
        }
        else {
          dm[i,j] <- dm[i-1,j]	
        }
      }
      # else if individual is removed stays removed
      else if (dm[i-1,j] == "R") {
        dm[i,j] <- dm[i-1,j]
      }
    }
  }
  return(as.data.frame(dm))
}


#------------------------------
# the stochastic SIR function
sir_stoch <- function(N, t, b, a) {
  df <- SIR(N, t, b, a)
  sir_m <- matrix(nrow = t, ncol = 3)
  colnames(sir_m) <- c("S","I","R")
  rownames(sir_m) <- c(1:t)
  
  sir_m[, "S"] <- apply(df == "S", 1, sum)
  sir_m[, "I"] <- apply(df == "I", 1, sum)
  sir_m[, "R"] <- apply(df == "R", 1, sum)
  
  return(as.data.frame(sir_m))
}


#------------------------------
# plot results of one SIR run
plot_sir <- function(df) {
  # takes input matrix with cols of iterations, S, I, R
  # and creates a plot with all three lines
  # S = black, I = red, R = blue
  t <- c(1:nrow(df))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(xpd=T, mar=par()$mar+c(0,0,0,3))
  matplot(t, df$S, main = 'SIR Cases', xlab='time steps', ylab='number of individuals', type='l', col=1, ylim = c(-1,max(df$S)))
  matlines(t, df$I, type = 'l', col = 'red')
  matlines(t, df$R, type = 'l', col = 'blue')
  legend(1.07*nrow(df), max(df$S), c('S', 'I', 'R'), cex = 0.8, col = c('black', 'red', 'blue'), lty=c(1, 1, 1))
  
  return(invisible(NULL)) #just to keep things honest
}


#------------------------------
# plot additional SIR runs, adding to an existing plot
add_plot_sir <- function(df) {
  # overlay additional SIR simulations onto an existing plot
  t <- c(1:nrow(df))
  matlines(t, df$S, col=1)
  matlines(t, df$I, col=2)
  matlines(t, df$R, col=3)
}


#------------------------------
# plot results of one stochastic SIR run
multi_plot_sir <- function(N, t, b, a, reps) {
  # run x number of stochastic SIR simulations and plot
  plot_sir(sir_stoch(N, t, b, a))	
  for(i in 2:reps) {
    add_plot_sir(sir_stoch(N, t, b, a))	
  }
  return(invisible(NULL)) #just to keep things honest
}


#------------------------------
# plot results of summarize_sir
error_plot_sir <- function(df){
  # takes the output from the super_summary and plots
  # the S,I,R classes with sd error bars for each timepoint
  
  t <- c(1:nrow(df))
  N <- sum(df[1,])
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(xpd=T, mar=par()$mar+c(0,0,0,3))
  
  matplot(t, df$muS,  main = 'SIR Cases', xlab='time steps', ylab='number of individuals', type='l', col=1, ylim = c(-1.5,max(df$muS)+2.5))
  segments(t, apply(df, 1, function(x) {if (df$muS-df$sdS < 0) {0} else {df$muS-df$sdS}}), t, apply(df, 1, function(x) {if (df$muS+df$sdS > N) {N} else {df$muS+df$sdS}}))
  
  matlines(t, df$muI, type = 'l', col = 'red')
  segments(t, df$muI-df$sdI, t, df$muI+df$sdI, col = 'red')
  
  matlines(t, df$muR, type = 'l', col = 'blue')
  segments(t, df$muR-df$sdR, t, df$muR+df$sdR, col = 'blue')
  
  legend(1.07*nrow(df), max(df$muS), c('S', 'I', 'R'), cex = 0.8, col = c('black', 'red', 'blue'), lty=c(1, 1, 1))
  
  return(invisible(NULL)) #just to keep things honest
}


#------------------------------
# the super_summary function:
summarize_sir <- function(N, t, b, a, reps) {
  # creates a series of stochasticSIR reps that
  # are summarized by the summarize_run function and
  # then loads the S,I,R counts from each timestep into 
  # separate S/I/R-specific matrices, and then the
  # S/I/R-specific matrices are used to create a 
  # super_summary matrix that has the mean and sd
  # of the counts of S/I/R for each timestep across
  # all reps
  
  # creates separate matrices to keep the SIR values for each timestep and each run
  S_val <- matrix(nrow=t, ncol=reps)
  I_val <- matrix(nrow=t, ncol=reps)
  R_val <- matrix(nrow=t, ncol=reps)
  
  # runs a stochastic sim, summarizes it, and loads the values for timesteps of that run into the summary matrices
  for(i in 1:reps) {
    df <- sir_stoch(N, t, b, a)
    #print(i)
    S_val[, i] <- df$S
    I_val[, i] <- df$I
    R_val[, i] <- df$R
  }
  
  # creates an overall summary matrix for the average and stdev for SIR values at each timestep across all runs
  df <- as.data.frame(matrix(nrow = t, ncol = 6))
  colnames(df) <- c("muS","sdS","muI","sdI","muR","sdR")
  rownames(df) <- c(1:t)
  
  # creates vectors that hold the means/sds for SIR values applying fn going down rows of the S/I/R value matrices
  S_rowmeans <- apply(S_val, 1, mean)
  S_rowstdev <- apply(S_val, 1, sd)
  I_rowmeans <- apply(I_val, 1, mean)
  I_rowstdev <- apply(I_val, 1, sd)
  R_rowmeans <- apply(R_val, 1, mean)
  R_rowstdev <- apply(R_val, 1, sd)
  
  # assign the values from the mean/sd vectors into the summary matrix
    df$muS <- round(S_rowmeans, 1)
    df$sdS <- round(S_rowstdev, 2)
    df$muI <- round(I_rowmeans, 1)
    df$sdI <- round(I_rowstdev, 2)
    df$muR <- round(R_rowmeans, 1)
    df$sdR <- round(R_rowstdev, 2)
    
  return(df)
}

