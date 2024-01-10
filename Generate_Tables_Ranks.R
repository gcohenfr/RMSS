#--------------------------
# GENERATE RANK TABLES 
#--------------------------

# Clear Memory
rm(list=ls())

# Required Library

# Required Source Files

#--------------------------

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

# Metric of interest
metric <- c("MSPE", "RC", "PR")[1]

# Signal-to-noise ratios considered
snr <- c(0.5, 1, 2)

# Correlation values 
rho_within <- c(0.8)

# Contamination proportions considered
tau_fixed <- c(0, 0.15, 0.3)[2]

# Sparsities considered
sparsity <- c(0.1, 0.2, 0.4)

# Ranking table
ranks <- matrix(0, nrow = 6, ncol = 9)
rownames(ranks) <- c("EN", "PENSE", "HuberEN", "SparseLTS", "RBSS", "RMSS")

# Current case ID
case_id <- 1

# Looping over the contamination scenarios
for(cont_scenario_val in cont_scenario){
  
  for(snr_val in snr){
    
    # Looping over the correlation values considered
    for(rho_within_val in rho_within){
      
      # Load the results
      filename <- paste0("results/results_n=", n, "_p=", p, 
                         "_cont=", cont_scenario_val, 
                         "_snr=", snr_val, "_rho=", rho_within_val,
                         "_tau=", tau_fixed, 
                         ".Rdata")
      load(filename)
      
      # Filling rank matrix
      for(sparsity_level in 1:length(sparsity)){
        
        rank_result <- rownames(output[[sparsity_level]])[
          order(apply(output[[sparsity_level]], 1:2, mean, na.rm = TRUE)[c(1:6), which(metric == c("MSPE", "RC", "PR"))],
                                                                decreasing = ifelse(metric != "MSPE", TRUE, FALSE))]
        ranks["EN", case_id] <- which("EN" == rank_result)
        ranks["PENSE", case_id] <- which("PENSE" == rank_result)
        ranks["HuberEN", case_id] <- which("Huber-EN" == rank_result)
        ranks["SparseLTS", case_id] <- which("SparseLTS" == rank_result)
        ranks["RBSS", case_id] <- which("RBSS" == rank_result)
        ranks["RMSS", case_id] <- which("RMSS" == rank_result)
        case_id <- case_id + 1
      }
    } 
  }
}

# Average rank
round(apply(ranks, 1, mean), 1)

# Lowest rank
round(apply(ranks, 1, max), 1)

