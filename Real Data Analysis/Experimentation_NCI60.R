#--------------------------------
# GENERATE OUTPUT - NCI-60 DATA
#--------------------------------

# Clear Memory
rm(list=ls())

# Required Library
library(parallel)

# Required Source Files
library(robustHD)

#--------------------------------

# Loading data
data("nci60")

# ______________
# Data analysis
# ______________

# Data dimension
dim(gene)
dim(protein)

# Protein MAD order
sort(apply(protein, 2, mad), decreasing = TRUE)

# Data pre-processing
y <- protein[, 159] 
correlations <- apply(gene, 2, corHuber, y)
keep <- partialOrder(abs(correlations), 500, decreasing = TRUE)
x <- scale(gene[, keep])
y <- scale(y)

# Outliers analysis
ddc_output <- cellWise::DDC(cbind(x, y), DDCpars = list(fastDDC = TRUE, silent = TRUE))
cellWise::cellMap(ddc_output$stdResid)

# Training data
train_id <- sample(1:nrow(x), 30, replace = FALSE)
x_train <- x[train_id,]
y_train <- y[train_id]

# Test data
x_test <- x[-train_id,]
y_test <- y[-train_id]

# EN
en <- tryCatch({
  
  cpu_en <- system.time(en_fit <- glmnet::cv.glmnet(x = x_train,
                                                    y = y_train,
                                                    alpha = 3/4))
  en_coef <- coef(en_fit, lambda = "lambda.min")[-1]
  preds_en <- predict(en_fit, x_test, lambda = "lambda.min")
  mspe_en <- robustbase::scaleTau2(abs(preds_en - y_test), mu.too = TRUE)[1]
  c(mspe_en, cpu_en["elapsed"])
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
  MSPE_sparseLTS <- robustbase::scaleTau2(abs(predict(sparseLTS_output, x_test) - y_test), mu.too = TRUE)
  CPU_sparseLTS <- sparseLTS_cpu["elapsed"]
  
  c(MSPE_sparseLTS, CPU_sparseLTS)
  
}, error = function(e){
  return(c(NA, NA))
}) 

# RMSS/RBSS
rbss_rmss <- tryCatch({
  
  # Number of models 
  n_models <- 5
  
  # h_grid
  h_grid <- ceiling(nrow(x_train) * 0.75)
  
  # RMSS
  cpu_rmss <- system.time(rmss_fit <- RMSS::cv.RMSS(x = x_train, y = y_train,
                                                    n_models = n_models,
                                                    h_grid = h_grid, t_grid = c(6, 8, 10), u_grid = c(1:n_models),
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
  mspe_rmss <- robustbase::scaleTau2(abs(y_test - preds_rmss), mu.too = TRUE)[1]
  
  # RBSS
  preds_rbss <- predict(rmss_fit, x_test, 
                        h_ind = rmss_fit$rbss_h_opt, 
                        t_ind = rmss_fit$rbss_t_opt, 
                        u_ind = rmss_fit$n_models)
  mspe_rbss <- robustbase::scaleTau2(abs(y_test - preds_rbss), mu.too = TRUE)[1]
  
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
  mspe_rf <- robustbase::scaleTau2(abs(y_test - preds_rf), mu.too = TRUE)[1]
  c(mspe_rf, cpu_rf["elapsed"])
}, error = function(e) {
  return(c(NA, NA))
})

en
sparselts
rbss_rmss
rf

sum(coef(sparseLTS_output)[-1] != 0)
sum(coef(rmss_fit)[-1] != 0)


















