#--------------------------------------
# GENERATE OUTPUT FROM THE SIMULATIONS
#--------------------------------------

# Clear Memory
rm(list=ls())

# Required Library
library(parallel)

# Required Source Files
source("Simulation_Run.R")

#--------------------------------------

# _______________________
# Setting the parameters
# _______________________

# Number of replications per simulation scenario
N <- 50
# Sample size for test data
m <- 2e3

# (*) Number of parameters and observations
p <- 500
n <- 50

# (*Adjust*) Covariance Scenarios
cont_scenario <- c("Cohen-Freue", "Thompson")

# Signal-to-noise ratios considered
snr <- c(0.5, 1, 2)

# Correlation values considered (within and between blocks)
group_size <- 25
rho_within <- c(0.8)
rho_between <- c(0.2)

# Contamination proportions considered
tau <- c(0, 0.1, 0.2, 0.3)

# Sparsities considered
sparsity <- c(0.1, 0.2, 0.4)

# Contamination parameters
k_lev <- 2
k_slo <- 100
cont_pred <- 0.2

# Number of models for RMSS
n_models <- 10

# Seed
seed <- 0

# Looping over the contamination scenarios
for(cont_scenario_val in cont_scenario){
  
  for(snr_val in snr){
    
    # Looping over the correlation values considered
    for(rho_within_val in rho_within){
      
      # Print simulation information
      cat("\n", "Contamination Scenario: ", cont_scenario_val)
      cat("\n", "SNR: ", snr_val)
      cat("\n", "Within-Block Correlation: ", rho_within_val, "\n")
      
      # Saving the file names for:
      # 1) The number of observations n
      # 2) The number of covariates p
      # 3) The contamination scenario cont_scenario
      # 4) The level of within-block correlation rho_within
      filename <- paste0("results/results_n=", n, "_p=", p, 
                         "_cont=", cont_scenario_val, 
                         "_snr=", snr_val, "_rho=", rho_within_val, ".Rdata")
      
      # Creating list for output
      output <- lapply(1:length(tau), function(t1) return(lapply(1:length(sparsity), function(t2) return(list()))))
      
      # Generate the output for the simulations
      for(tau_ind in 1:length(tau)){
        for(sparsity_ind in 1:length(sparsity)){
          
          # Print tau and sparsity
          cat("\n", "tau: ", tau[tau_ind])
          cat("\n", "sparsity: ", sparsity[sparsity_ind], "\n")
          
          # A three-dimensional array
          # MSPE, RC, PR, CPU for each method over the replications
          output[[tau_ind]][[sparsity_ind]] <- Simulation_Run(N = N, m = m, n = n, p = p, 
                                                              cont_scenario = cont_scenario_val, 
                                                              snr = snr_val, 
                                                              rho_within = rho_within_val, rho_between = rho_between, 
                                                              group_size = group_size, 
                                                              tau = tau[tau_ind], sparsity = sparsity[sparsity_ind], 
                                                              k_lev = k_lev, k_slo = k_slo, cont_pred = cont_pred,
                                                              n_models,
                                                              seed = seed)
        }
      } 

      # Saving the file on the git repository
      save.image(filename)
    }
  }
}


