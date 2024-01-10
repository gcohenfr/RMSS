#----------------------------
# GENERATE OUTPUT - BBS DATA
#----------------------------

# Clear Memory
rm(list=ls())

# Required Library
library(parallel)

# Required Source Files
library(abess)

#----------------------------

# ______________
# Loading Data
# ______________

# Loading data
data("trim32")

# Processsing data
y <- scale(trim32[, 1])
x <- scale(trim32[, -1])

# ____________________________
# Fixed Simulation Parameters
# ____________________________

# Seed for simulation
seed <- 1
# Number of replications
N <- 50

# _______________________________
# Variable Simulation Parameters
# _______________________________

# Setting the sample size of the training data
n <- c(50)

# Contamination proportion
tau <- 0.25
# Contamination of predictors proportion
cont_pred <- 0.2
n_cont <- floor(cont_pred*ncol(x))
# Contamination cell mean
cont_mean <- 25

# Contamination of response
contaminate_y <- c(TRUE)

# ___________________
# Running Simulation
# ___________________

for(n_size in n){
  
  # Print current sample size
  cat("n: ", n_size, "\n")
  
  for(y_cont in contaminate_y){
    
    # Print current contamination status
    cat("contaminate_y: ", y_cont, "\n")
    
    # Array to store prediction output
    predictions_output <- array(dim = c(9, 3, N))
    rownames(predictions_output) <- c("EN", "PENSE", "HuberEN", "SparseLTS", 
                                      "RBSS", "RMSS", "RMSS_DIV",
                                      "RF", "RGLM")
    colnames(predictions_output) <- c("MSPE", "AMSPE", "CPU")
    
    # List to store variable selections
    selections <- lapply(1:9, function(t) list())
    names(selections) <- c("EN", "PENSE", "HuberEN", "SparseLTS", 
                           "RBSS", "RMSS", "RMSS_DIV",
                           "RF", "RGLM")
    
    # Setting the seed 
    set.seed(seed)
    
    # Generate training IDs
    train_id <- lapply(1:N, function(t, x) sample(1:nrow(x), n_size, replace = FALSE), x = x)
    
    for(iter in 1:N){
      
      # Print current iteration
      cat("iter: ", iter, "\n")
      
      # Training data
      x_train <- as.matrix(x[train_id[[iter]],])
      y_train <- y[train_id[[iter]]]
      
      # Contamination of training data
      cont_ind <- sample(1:n_size, floor(tau*n_size))
      for(cont_id in cont_ind)
        x_train[cont_id, sample(1:ncol(x_train), n_cont)] <- rnorm(n_cont, mean = cont_mean)
      if(y_cont)
        y_train[cont_ind] <- rnorm(length(cont_ind), mean = cont_mean)
      
      # Test data
      x_test <- as.matrix(x[-train_id[[iter]],])
      y_test <- y[-train_id[[iter]]]
      
      # EN
      en <- tryCatch({
        
        cpu_en <- system.time(en_fit <- glmnet::cv.glmnet(x = x_train,
                                                          y = y_train,
                                                          alpha = 3/4))
        preds_en <- predict(en_fit, x_test, lambda = "lambda.min")
        mspe_en <- mean((preds_en - y_test)^2)
        
        selections$EN[[iter]] <- which(coef(en_fit, s = "lambda.min")[-1] != 0)
        
        c(mspe_en, NA, cpu_en["elapsed"])
      }, error = function(e) {
        selections$EN[[iter]] <- NA
        return(c(NA, NA, NA))
      })
      predictions_output["EN",, iter] <- en
      
      # PENSE
      pense <- tryCatch({
        
        cpu_pense <- system.time(pense_fit <- pense::adapense_cv(x_train, y_train,
                                                                 alpha = 3/4, cv_k = 5, cv_repl = 1, cl = NULL,
                                                                 eps = 1e0, explore_tol = 1e3,
                                                                 enpy_opts = pense::enpy_options(retain_max = 5)))
        preds_pense <- predict(pense_fit, newdata = x_test, lambda = pense_fit$lambda[[1]][which.min(pense_fit$cvres$cvavg)])
        mspe_pense <- mean((y_test - preds_pense)^2)
        
        selections$PENSE[[iter]] <- as.numeric(which(coef(pense_fit, lambda = pense_fit$lambda[[1]][which.min(pense_fit$cvres$cvavg)])[-1] != 0))
        
        c(mspe_pense, NA, cpu_pense["elapsed"])
      }, error = function(e) {
        selections$PENSE[[iter]] <- NA
        return(c(NA, NA, NA))
      })
      predictions_output["PENSE",, iter] <- pense
      
      # Huber-EN
      huber <- tryCatch({
        
        cpu_huber <- system.time(huber_fit <- hqreg::cv.hqreg(x_train, y_train, 
                                                              method = "huber", alpha = 3/4,
                                                              nlambda = 50,
                                                              nfolds = 5))
        preds_huber <- predict(huber_fit, x_test)
        mspe_huber <- mean((y_test - preds_huber)^2)
        
        selections$HuberEN[[iter]] <- as.numeric(which(coef(huber_fit)[-1] != 0))
        
        c(mspe_huber, NA, cpu_huber["elapsed"])
      }, error = function(e) {
        selections$HuberEN[[iter]] <- NA
        return(c(NA, NA, NA))
      })
      predictions_output["HuberEN",, iter] <- huber
      
      # Sparse LTS
      sparselts <- tryCatch({
        
        lambda_max <- robustHD::lambda0(x_train, y_train)
        lambda_grid <- rev(exp(seq(log(1e-2*lambda_max), log(lambda_max), length = 50)))
        sparseLTS_cpu <- system.time(
          sparseLTS_output <- robustHD::sparseLTS(x = x_train, 
                                                  y = y_train, 
                                                  lambda = lambda_grid,
                                                  mode = "lambda",
                                                  tol = 1e-2,
                                                  cluster = NULL,
                                                  ncores = 5))
        MSPE_sparseLTS <- mean((predict(sparseLTS_output, x_test) - y_test)^2)
        CPU_sparseLTS <- sparseLTS_cpu["elapsed"]
        
        selections$SparseLTS[[iter]] <- as.numeric(which(coef(sparseLTS_output)[-1] != 0))
        
        c(MSPE_sparseLTS, NA, CPU_sparseLTS)
      }, error = function(e){
        selections$SparseLTS[[iter]] <- NA
        return(c(NA, NA, NA))
      }) 
      predictions_output["SparseLTS",, iter] <- sparselts
      
      # RMSS/RBSS
      rbss_rmss <- tryCatch({
        
        # Number of models 
        n_models <- 10
        
        # h_grid
        h_grid <- ceiling(nrow(x_train) * 0.75)
        
        # RMSS
        cpu_rmss <- system.time(rmss_fit <- RMSS::cv.RMSS(x = x_train, y = y_train,
                                                          n_models = n_models,
                                                          h_grid = h_grid, 
                                                          t_grid = c(floor(0.3*nrow(x_train)), 
                                                                     floor(0.4*nrow(x_train)), 
                                                                     floor(0.5*nrow(x_train))), 
                                                          u_grid = c(1:n_models),
                                                          initial_estimator = "robStepSplitReg",
                                                          tolerance = 1e-1,
                                                          max_iter = 1e3,
                                                          neighborhood_search = FALSE,
                                                          neighborhood_search_tolerance = 1e-1,
                                                          cv_criterion = "tau",
                                                          n_folds = 5,
                                                          gamma = 1, 
                                                          n_threads = 5))
        preds_rmss <- predict(rmss_fit, x_test)
        mspe_rmss <- mean((y_test - preds_rmss)^2)
        
        preds_rmss_div <- predict(rmss_fit, x_test, t_ind = which.min(rmss_fit$cv_error[[1]]), u_ind = 1)
        mspe_rmss_div <- mean((y_test - preds_rmss_div)^2)
        
        mspe_models <- numeric(n_models)
        for(model_id in 1:n_models){
          
          mspe_models[model_id] <- mean((y_test - predict(rmss_fit, x_test, group_index = model_id))^2)
        }
        
        mspe_models_div <- numeric(n_models)
        for(model_id in 1:n_models){
          
          mspe_models_div[model_id] <- mean((y_test - predict(rmss_fit, x_test, group_index = model_id,
                                                              t_ind = which.min(rmss_fit$cv_error[[1]]), u_ind = 1))^2)
        }
        
        selections$RMSS[[iter]] <- table(apply(coef(rmss_fit, individual_models = TRUE)[-1,], 2, function(t) which(t != 0)))
        selections$RMSS_DIV[[iter]] <- as.numeric(apply(coef(rmss_fit, individual_models = TRUE,
                                                             t_ind = which.min(rmss_fit$cv_error[[1]]), u_ind = 1)[-1,], 2, 
                                                        function(t) which(t != 0)))
        
        # RBSS
        preds_rbss <- predict(rmss_fit, x_test, 
                              h_ind = rmss_fit$rbss_h_opt, 
                              t_ind = rmss_fit$rbss_t_opt, 
                              u_ind = rmss_fit$n_models)
        mspe_rbss <- mean((y_test - preds_rbss)^2)
        
        # Combined output
        rbind(c(mspe_rbss, NA, NA),
              c(mspe_rmss, mean(mspe_models), cpu_rmss["elapsed"]),
              c(mspe_rmss_div, mean(mspe_models_div), NA))
      }, error = function(e) {
        selections$RBSS[[iter]] <- NA
        selections$RMSS[[iter]] <- NA
        selections$RMSS_DIV[[iter]] <- NA
        return(rbind(c(NA, NA, NA),
                     c(NA, NA, NA),
                     c(NA, NA, NA)))
      })
      predictions_output["RBSS",, iter] <- rbss_rmss[1,]
      predictions_output["RMSS",, iter] <- rbss_rmss[2,]
      predictions_output["RMSS_DIV",, iter] <- rbss_rmss[3,]
      
      # Random Forest
      rf <- tryCatch({
        
        cpu_rf <- system.time(rf_fit <- randomForest::randomForest(x = x_train, 
                                                                   y = c(y_train)))
        preds_rf <- predict(rf_fit, newdata = x_test)
        mspe_rf <- mean((y_test - preds_rf)^2)
        
        mspe_trees <- numeric(n_models)
        mspe_trees[1] <- mean((y_test - predict(rf_fit, newdata = x_test))^2)
        for(tree_id in 2:rf_fit$ntree){
          
          rf_tree <- randomForest::randomForest(x = x_train, y = y_train, ntree = 1)
          mspe_trees[tree_id] <- mean((y_test - predict(rf_tree, newdata = x_test))^2)
          rf_fit <- randomForest::combine(rf_fit, rf_tree)
        }
        
        c(mspe_rf, mean(mspe_trees), cpu_rf["elapsed"])
      }, error = function(e) {
        return(c(NA, NA, NA))
      })
      predictions_output["RF",, iter] <- rf
      
      # Random GLM
      rglm <- tryCatch({
        
        cpu_rglm <- system.time(rglm_fit <- randomGLM::randomGLM(x = x_train, 
                                                                 y = y_train, 
                                                                 classify = FALSE,
                                                                 randomSeed = NULL))
        preds_rglm <- predict(rglm_fit, newdata = x_test)
        mspe_rglm <- mean((y_test - preds_rglm)^2)
        
        mspe_glm <- numeric(length(rglm_fit$models))
        for(model_id in 1:length(rglm_fit$models)){
          
          model_pred <- rglm_fit$interceptOfForwardRegression[[model_id]] + 
            x_test[, rglm_fit$featuresInForwardRegression[[model_id]]] %*% rglm_fit$coefOfForwardRegression[[model_id]]
          mspe_glm[model_id] <- mean((model_pred - y_test)^2)
        }
        
        c(mspe_rglm, mean(mspe_glm), cpu_rglm["elapsed"])
      }, error = function(e) {
        return(c(NA, NA, NA))
      })
      predictions_output["RGLM",, iter] <- rglm
    }
    
    filename <- paste0("results/results_BBS_n=", n_size, "_y_cont=", y_cont, ".RData")
    save.image(filename)
  }
}




















