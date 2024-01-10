#------------------------------
# GENERATE COMPUTING TABLE
#------------------------------

# Clear Memory
rm(list=ls())

# Required Library

# Required Source Files

#------------------------------

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
cont_scenario <- c("Cohen-Freue", "Thompson")[1]

# Signal-to-noise ratios considered
snr <- c(0.5, 1, 2)

# Correlation values considered (within and between blocks)
group_size <- 25
rho_within <- c(0.8)
rho_between <- c(0.2)

# Contamination proportions considered
tau <- c(0, 0.15, 0.3)

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

# Vector to store CPU times
cpu_time <- numeric(6)
names(cpu_time) <- c("EN", "PENSE", "HuberEN", "SparseLTS", "RBSS", "RMSS")

# Looping over the contamination scenarios
for(cont_scenario_val in cont_scenario){
  
  for(snr_val in snr){
    
    # Looping over the correlation values considered
    for(rho_within_val in rho_within){
      
      # Generate the output for the simulations
      for(tau_val in tau){
        

        # Load the results
        filename <- paste0("results/results_n=", n, "_p=", p, 
                           "_cont=", cont_scenario_val, 
                           "_snr=", snr_val, "_rho=", rho_within_val,
                           "_tau=", tau_val, 
                           ".Rdata")
        load(filename)
        
        # Adding CPU times
        for(sparsity_ind in 1:length(sparsity)){
          
          cpu_time <- cpu_time + apply(output[[sparsity_ind]], 1:2, mean)[1:6, "CPU"] * N
        }
      } 
    }
  }
}

# CPU time 
round(cpu_time / (N * 3^3), 1)

