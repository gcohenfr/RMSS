#---------------------------------------------
# DATA SIMULATION
#---------------------------------------------

# Function to simulate training and test data for all parameters specified
Data_Simulation <- function(N, m, n, p, cont_scenario, 
                            snr, rho_within, rho_between, group_size, 
                            tau, sparsity, 
                            k_lev, k_slo, cont_pred,
                            seed) {
  
  # Setting the seed 
  set.seed(0)
  
  # Number of nonzero regression coefficients
  p_active <- floor(p*sparsity)
  
  # Block Correlation
  sigma.mat <- matrix(0, p, p)
  sigma.mat[1:p_active, 1:p_active] <- rho_between
  for(group in 0:(p_active/group_size - 1))
    sigma.mat[(group*group_size+1):(group*group_size+group_size),(group*group_size+1):(group*group_size+group_size)] <- rho_within
  diag(sigma.mat) <- 1
  
  # Simulation of beta vector
  true.beta <- c(runif(p_active, 0, 5)*(-1)^rbinom(p_active, 1, 0.7), rep(0, p - p_active))
  
  # Setting the SD of the variance
  sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))
  
  # List to store clean data
  clean_data <- list()
  
  # Simulation of uncontaminated data
  for(rep in 1:N){
    
    # Simulation of uncontaminated data 
    clean_data[[rep]] <- list()
    clean_data[[rep]]$x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
    clean_data[[rep]]$y <- clean_data[[rep]]$x %*% true.beta + rnorm(n, 0, sigma)
  }
  
  # List to store training and test data
  train_data <- list()
  
  # Data contamination
  if(tau == 0){
    
    train_data <- clean_data
  } else{
    
    # Simulation of training data
    if(cont_scenario == "Cohen-Freue"){
      
      # Setting contaminated coefficients
      beta_cont <- true.beta
      beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
      beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
      
      for(rep in 1:N){
        
        # Contamination of data 
        x_train <- clean_data[[rep]]$x
        y_train <- clean_data[[rep]]$y
        contamination_indices <- 1:floor(n*tau)
        for(cont_id in contamination_indices){
          
          a <- runif(p, min = -1, max = 1)
          a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
          x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
          y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
        }
        
        # Store training data
        train_data[[rep]] <- list()
        train_data[[rep]]$x <- x_train
        train_data[[rep]]$y <- y_train
      }
      
    } else if (cont_scenario == "Thompson"){
      
      # Proportion of contaminated predictors for contaminated samples
      n_cont <- floor(cont_pred*p)
      
      for(rep in 1:N){
        
        # Contamination of data 
        contamination_indices <- 1:floor(n*tau)
        x_train <- clean_data[[rep]]$x
        y_train <- clean_data[[rep]]$y
        for(cont_id in contamination_indices)
          x_train[cont_id, sample(1:ncol(x_train), n_cont)] <- rnorm(n_cont, mean = 10)
        y_train[contamination_indices] <- rnorm(length(contamination_indices), mean = 10*sigma, sd = sigma)
        
        # Store training data
        train_data[[rep]] <- list()
        train_data[[rep]]$x <- x_train
        train_data[[rep]]$y <- y_train
      }
    }
  }
  
  # Simulation of test data
  x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
  y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)
  
  # List with simulated data and parameters
  sim_data <- list(N = N, m = m, n = n, p = p,
                   train_data = train_data, 
                   clean_data = clean_data,
                   test_data = list(x = x_test, y = y_test),
                   sigma = sigma,
                   active_predictors = 1:p_active,
                   tau = tau)
  return(sim_data)
}




