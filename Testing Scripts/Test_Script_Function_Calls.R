#--------------------------------------
# TEST SCRIPT - FUNCTION CALLS
#--------------------------------------

# Clear Memory
rm(list=ls())

# Required library
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
cont_scenario <- c("Cohen-Freue", "Thompson")[1]

# Signal-to-noise ratios considered
snr <- c(0.5, 1, 2)[1]

# Correlation values considered (within and between blocks)
group_size <- 25
rho_within <- c(0.8)
rho_between <- c(0.2)

# Contamination proportions considered
tau <- c(0, 0.1, 0.2, 0.3)[2]

# Sparsities considered
sparsity <- c(0.1, 0.2, 0.4)[2]

# Contamination parameters
k_lev <- 2
k_slo <- 100
cont_pred <- 0.2

# Number of models for RMSS
n_models <- 10

# Seed
seed <- 0

# ____________________
# Simulating the data
# ____________________

# Simulation of training and test data
simulated_data <- Data_Simulation(N, m, n, p, cont_scenario, 
                                  snr, rho_within, rho_between, group_size, tau, sparsity, 
                                  k_lev, k_slo, cont_pred,
                                  seed)

# ___________________________
# Generating the predictions
# ___________________________

# Array to store output
predictions_output <- array(dim = c(8, 4, simulated_data$N))
rownames(predictions_output) <- c("EN", "PENSE", "Huber-EN", "SparseLTS", 
                                  "RBSS", "RMSS", "RF", "RGLM")
colnames(predictions_output) <- c("MSPE", "RC", "PR", "CPU")

# Setting trimming grid for RMSS
h_grid <- round(c(1 - simulated_data$tau - 0.1, 1 - simulated_data$tau - 0.05, 1 - simulated_data$tau) * simulated_data$n)

# Cluster for PENSE and SparseLTS
cluster <- makeCluster(5)

# Looping over training sets
for(rep in 1:simulated_data$N){
  
  # Print iteration
  cat("\n", "Iteration: ", rep)
  
  # EN
  cat("\n", "Computing EN.")
  en <- tryCatch({
    
    cpu_en <- system.time(en_fit <- glmnet::cv.glmnet(x = simulated_data$train_data[[rep]]$x,
                                                      y = simulated_data$train_data[[rep]]$y,
                                                      alpha = 3/4))
    en_coef <- coef(en_fit, lambda = "lambda.min")[-1]
    recall_en <- sum(which((en_coef[-1]!=0)) %in% simulated_data$active_predictors)/length(simulated_data$active_predictors)
    if(sum(en_coef[-1]!=0) > 0)
      prec_en <- sum(which((en_coef[-1]!=0)) %in% simulated_data$active_predictors)/sum(en_coef[-1]!=0) else
        prec_en <- 0
    preds_en <- predict(en_fit, simulated_data$test_data$x, lambda = "lambda.min")
    mspe_en <- mean((preds_en - simulated_data$test_data$y)^2)/simulated_data$sigma^2
    c(mspe_en, recall_en, prec_en, cpu_en["elapsed"])
  }, error = function(e) {
    return(c(NA, NA, NA, NA))
  })
  predictions_output["EN",, rep] <- en
  
  # PENSE
  cat("\n", "Computing PENSE.")
  pense <- tryCatch({
    
    cpu_pense <- system.time(pense_fit <- pense::adapense_cv(simulated_data$train_data[[rep]]$x, simulated_data$train_data[[rep]]$y,
                                                             alpha = 3/4, cv_k = 5, cv_repl = 1, cl = cluster,
                                                             eps = 1e0, explore_tol = 1e3,
                                                             enpy_opts = pense::enpy_options(retain_max = 5)))
    pense_coef <- coef(pense_fit, lambda = pense_fit$lambda[[1]][which.min(pense_fit$cvres$cvavg)])
    recall_pense <- sum(which((pense_coef[-1]!=0)) %in% simulated_data$active_predictors)/length(simulated_data$active_predictors)
    if(sum(pense_coef[-1]!=0) > 0)
      prec_pense <- sum(which((pense_coef[-1]!=0)) %in% simulated_data$active_predictors)/sum(pense_coef[-1]!=0) else
        prec_pense <- 0
    preds_pense <- predict(pense_fit, newdata = simulated_data$test_data$x, lambda = pense_fit$lambda[[1]][which.min(pense_fit$cvres$cvavg)])
    mspe_pense <- mean((simulated_data$test_data$y - preds_pense)^2)/simulated_data$sigma^2
    c(mspe_pense, recall_pense, prec_pense, cpu_pense["elapsed"])
  }, error = function(e) {
    return(c(NA, NA, NA, NA))
  })
  predictions_output["PENSE",, rep] <- pense
  
  # Huber-EN
  cat("\n", "Computing Huber-EN.\n")
  huber <- tryCatch({
    
    cpu_huber <- system.time(huber_fit <- hqreg::cv.hqreg(simulated_data$train_data[[rep]]$x, simulated_data$train_data[[rep]]$y, 
                                                          method = "huber", alpha = 3/4,
                                                          nlambda = 50,
                                                          nfolds = 5))
    huber_coef <- coef(huber_fit)
    recall_huber <- sum(which((huber_coef[-1]!=0))  %in% simulated_data$active_predictors)/length(simulated_data$active_predictors)
    if(sum(huber_coef[-1]!=0) > 0)
      prec_huber <- sum(which((huber_coef[-1]!=0))  %in% simulated_data$active_predictors)/sum(huber_coef[-1]!=0) else
        prec_huber <- 0
    preds_huber <- predict(huber_fit, simulated_data$test_data$x)
    mspe_huber <- mean((simulated_data$test_data$y - preds_huber)^2)/simulated_data$sigma^2
    c(mspe_huber, recall_huber, prec_huber, cpu_huber["elapsed"])
  }, error = function(e) {
    return(c(NA, NA, NA, NA))
  })
  predictions_output["Huber-EN",, rep] <- huber
  
  # Sparse LTS
  cat("\n", "Computing SparseLTS.")
  sparselts <- tryCatch({
    
    lambda_max <- robustHD::lambda0(simulated_data$train_data[[rep]]$x, simulated_data$train_data[[rep]]$y)
    lambda_grid <- rev(exp(seq(log(1e-2*lambda_max), log(lambda_max), length = 50)))
    cpu_sparselts <- system.time(sparselts_fit <- robustHD::sparseLTS(simulated_data$train_data[[rep]]$x, c(simulated_data$train_data[[rep]]$y),
                                                                      mode = "lambda", lambda = lambda_grid, tol = 1e-1,
                                                                      cluster = cluster,
                                                                      ncores = 5))
    sparselts_coef <- as.numeric(coef(sparselts_fit))
    recall_sparselts <- sum(which((sparselts_coef[-1]!=0)) %in% simulated_data$active_predictors)/length(simulated_data$active_predictors)
    if(sum(sparselts_coef[-1]!=0) > 0)
      prec_sparselts <- sum(which((sparselts_coef[-1]!=0)) %in% simulated_data$active_predictors)/sum(sparselts_coef[-1]!=0) else
        prec_sparselts <- 0
    preds_sparselts <- predict(sparselts_fit, simulated_data$test_data$x)
    mspe_sparselts <- mean((simulated_data$test_data$y - preds_sparselts)^2)/simulated_data$sigma^2
    c(mspe_sparselts, recall_sparselts, prec_sparselts, cpu_sparselts["elapsed"])
  }, error = function(e) {
    return(c(NA, NA, NA, NA))
  })
  predictions_output["SparseLTS",, rep] <- sparselts
  
  # RMSS/RBSS
  cat("\n", "Computing RMSS/RBSS")
  rbss_rmss <- tryCatch({
    
    # RMSS
    cpu_rmss <- system.time(rmss_fit <- RMSS::cv.RMSS(x = simulated_data$train_data[[rep]]$x, y = simulated_data$train_data[[rep]]$y,
                                                      n_models = n_models,
                                                      h_grid = h_grid, t_grid = c(10, 15, 20), u_grid = c(1:n_models),
                                                      initial_estimator = "robStepSplitReg",
                                                      tolerance = 1e-1,
                                                      max_iter = 1e3,
                                                      neighborhood_search = FALSE,
                                                      neighborhood_search_tolerance = 1e-1,
                                                      cv_criterion = "tau",
                                                      n_folds = 5,
                                                      gamma = 1, 
                                                      n_threads = 5))
    rmss_coef <- as.numeric(coef(rmss_fit))
    recall_rmss <- sum(which((rmss_coef[-1]!=0)) %in% simulated_data$active_predictors)/length(simulated_data$active_predictors)
    if(sum(rmss_coef[-1]!=0) > 0)
      prec_rmss <- sum(which((rmss_coef[-1]!=0)) %in% simulated_data$active_predictors)/sum(rmss_coef[-1]!=0) else
        prec_rmss <- 0
    preds_rmss <- predict(rmss_fit, simulated_data$test_data$x)
    mspe_rmss <- mean((simulated_data$test_data$y - preds_rmss)^2)/simulated_data$sigma^2
    
    # RBSS
    rbss_coef <- as.numeric(coef(rmss_fit, 
                                 h_ind = rmss_fit$rbss_h_opt, 
                                 t_ind = rmss_fit$rbss_t_opt, 
                                 u_ind = rmss_fit$n_models))
    recall_rbss <- sum(which((rbss_coef[-1]!=0)) %in% simulated_data$active_predictors)/length(simulated_data$active_predictors)
    if(sum(rbss_coef[-1]!=0) > 0)
      prec_rbss <- sum(which((rbss_coef[-1]!=0)) %in% simulated_data$active_predictors)/sum(rbss_coef[-1]!=0) else
        prec_rbss <- 0
    preds_rbss <- predict(rmss_fit, simulated_data$test_data$x, 
                          h_ind = rmss_fit$rbss_h_opt, 
                          t_ind = rmss_fit$rbss_t_opt, 
                          u_ind = rmss_fit$n_models)
    mspe_rbss <- mean((simulated_data$test_data$y - preds_rbss)^2)/simulated_data$sigma^2
    
    # Combined output
    rbind(c(mspe_rbss, recall_rbss, prec_rbss, NA),
          c(mspe_rmss, recall_rmss, prec_rmss, cpu_rmss["elapsed"]))
  }, error = function(e) {
    return(rbind(c(NA, NA, NA, NA),
                 c(NA, NA, NA, NA)))
  })
  predictions_output["RBSS",, rep] <- rbss_rmss[1,]
  predictions_output["RMSS",, rep] <- rbss_rmss[2,]
  
  # Random Forest
  cat("\n", "Computing RF.")
  rf <- tryCatch({
    
    cpu_rf <- system.time(rf_fit <- randomForest::randomForest(x = simulated_data$train_data[[rep]]$x, 
                                                               y = c(simulated_data$train_data[[rep]]$y)))
    preds_rf <- predict(rf_fit, newdata = simulated_data$test_data$x)
    mspe_rf <- mean((simulated_data$test_data$y - preds_rf)^2)/simulated_data$sigma^2
    c(mspe_rf, NA, NA, cpu_rf["elapsed"])
  }, error = function(e) {
    return(c(NA, NA, NA, NA))
  })
  predictions_output["RF",, rep] <- rf
  
  # Random GLM
  cat("\n", "Computing RGLM.\n")
  rglm <- tryCatch({
    
    cpu_rglm <- system.time(rglm_fit <- randomGLM::randomGLM(x = simulated_data$train_data[[rep]]$x, 
                                                             y = simulated_data$train_data[[rep]]$y, 
                                                             classify = FALSE, 
                                                             nThreads = 5))
    preds_rglm <- predict(rglm_fit, newdata = simulated_data$test_data$x)
    mspe_rglm <- mean((simulated_data$test_data$y - preds_rglm)^2)/simulated_data$sigma^2
    c(mspe_rglm, NA, NA, cpu_rglm["elapsed"])
  }, error = function(e) {
    return(c(NA, NA, NA, NA))
  })
  predictions_output["RGLM",, rep] <- rglm
}

# Ending the cluster
stopCluster(cluster)




















