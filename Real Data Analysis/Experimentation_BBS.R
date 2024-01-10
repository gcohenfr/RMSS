#--------------------------------
# GENERATE OUTPUT - NCI-60 DATA
#--------------------------------

# Clear Memory
rm(list=ls())

# Required Library
library(parallel)

# Required Source Files
library(abess)

#--------------------------------

# Loading data
data("trim32")

# ______________
# Data analysis
# ______________

# Processsing data
y <- scale(trim32[, 1])
x <- scale(trim32[, -1])

# Outliers analysis
ddc_output <- cellWise::DDC(cbind(x, y), DDCpars = list(fastDDC = TRUE, silent = TRUE))
cellWise::cellMap(ddc_output$stdResid)

# Setting the sample size of the training data
n <- 50

# Contamination proportion
tau <- 0.25
# Contamination of predictors proportion
cont_pred <- 0.2
n_cont <- floor(cont_pred*ncol(x))

# Training data
train_id <- sample(1:nrow(x), n, replace = FALSE)
x_train <- as.matrix(x[train_id,])
y_train <- y[train_id]

# Contamination of training data
cont_ind <- sample(1:n, floor(tau*n))
for(cont_id in cont_ind)
  x_train[cont_id, sample(1:ncol(x_train), n_cont)] <- rnorm(n_cont, mean = 15)
# y_train[cont_ind] <- rnorm(length(cont_ind), mean = 15)

# Test data
x_test <- as.matrix(x[-train_id,])
y_test <- y[-train_id]

# EN
en <- tryCatch({
  
  cpu_en <- system.time(en_fit <- glmnet::cv.glmnet(x = x_train,
                                                    y = y_train,
                                                    alpha = 3/4))
  en_coef <- coef(en_fit, lambda = "lambda.min")[-1]
  preds_en <- predict(en_fit, x_test, lambda = "lambda.min")
  mspe_en <- mean((preds_en - y_test)^2)
  c(mspe_en, cpu_en["elapsed"])
}, error = function(e) {
  return(c(NA, NA))
})

# PENSE
pense <- tryCatch({
  
  cpu_pense <- system.time(pense_fit <- pense::adapense_cv(x_train, y_train,
                                                           alpha = 3/4, cv_k = 5, cv_repl = 1, cl = NULL,
                                                           eps = 1e0, explore_tol = 1e3,
                                                           enpy_opts = pense::enpy_options(retain_max = 5)))
  preds_pense <- predict(pense_fit, newdata = x_test, lambda = pense_fit$lambda[[1]][which.min(pense_fit$cvres$cvavg)])
  mspe_pense <- mean((y_test - preds_pense)^2)
  
  c(mspe_pense, cpu_pense["elapsed"])
}, error = function(e) {
  return(c(NA, NA))
})

# Sparse LTS
sparselts <- tryCatch({
  
  lambda_max <- robustHD::lambda0(x_train, y_train)
  lambda_grid <- rev(exp(seq(log(1e-2*lambda_max), log(lambda_max), length = 10)))
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
  
  c(MSPE_sparseLTS, CPU_sparseLTS)
  
}, error = function(e){
  return(c(NA, NA))
}) 

# MSS
mss <- tryCatch({

  # Number of models
  n_models <- 10

  # MSS
  cpu_mss <- system.time(mss_fit <- PSGD::cv.PSGD(x = x_train, y = y_train, n_models = n_models,
                                                  model_type = c("Linear", "Logistic")[1], include_intercept = TRUE,
                                                  split = c(1, 2, 3, 4, 5),
                                                  size = c(floor(0.3*nrow(x_train)),
                                                           floor(0.4*nrow(x_train)),
                                                           floor(0.5*nrow(x_train))),
                                                  max_iter = 20,
                                                  cycling_iter = 0,
                                                  n_folds = 5))
  preds_psgd <- predict(mss_fit, newx = x_test, group_index = 1:n_models)
  mspe_psgd <- mean((y_test - preds_psgd)^2)
  cpu_mss <- cpu_mss["elapsed"]

  c(mspe_psgd, cpu_mss)
}, error = function(e) {
  return(c(NA))
})

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
  preds_rmss <- predict(rmss_fit, x_test, u_ind = 1)
  mspe_rmss <- mean((y_test - preds_rmss)^2)
  
  # RBSS
  preds_rbss <- predict(rmss_fit, x_test, 
                        h_ind = rmss_fit$rbss_h_opt, 
                        t_ind = rmss_fit$rbss_t_opt, 
                        u_ind = rmss_fit$n_models)
  mspe_rbss <- mean((y_test - preds_rbss)^2)
  
  # Combined output
  rbind(c(mspe_rbss, NA),
        c(mspe_rmss, cpu_rmss["elapsed"]))
}, error = function(e) {
  return(rbind(c(NA, NA),
               c(NA, NA)))
})

# Random Forest
rf <- tryCatch({
  
  cpu_rf <- system.time(rf_fit <- randomForest::randomForest(x = x_train, 
                                                             y = c(y_train)))
  preds_rf <- predict(rf_fit, newdata = x_test)
  mspe_rf <- mean((y_test - preds_rf)^2)
  c(mspe_rf, cpu_rf["elapsed"])
}, error = function(e) {
  return(c(NA, NA))
})

en
pense
sparselts
mss
rbss_rmss
rf

rmss_fit$u_opt
rmss_fit$t_opt

sum(coef(sparseLTS_output)[-1] != 0)
sum(coef(rmss_fit)[-1] != 0)


















