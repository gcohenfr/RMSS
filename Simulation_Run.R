#----------------------------
# RUN THE SIMULATIONS 
#----------------------------

# Required Source Files
source("Data_Simulation.R")
source("Generate_Predictions.R")

#----------------------------

# Function to run simulation for specified parameters
Simulation_Run <- function(N, m, n, p, cont_scenario, 
                           snr, rho_within, rho_between, group_size, tau, sparsity, 
                           k_lev, k_slo, cont_pred,
                           n_models,
                           seed){
  
  # Simulation of training and test data
  simulated_data <- Data_Simulation(N, m, n, p, cont_scenario, 
                                    snr, rho_within, rho_between, group_size, tau, sparsity, 
                                    k_lev, k_slo, cont_pred,
                                    seed)
  
  # Metrics for methods
  output <- Generate_Predictions(simulated_data,
                                 n_models)
  
  # Return the output
  return(output)
}
