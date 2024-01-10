######################################################################
##
# TRAINING: BiT1-discovery; TEST: BiT1-BiT2 
#   patients in the training are not included in the test set 
#   dataset built with "matched_training_test_BiT1_2.R" file
##
######################################################################
## Based on PENSE-PENSEM
#  Alpha=0.6 and lambda are selected based on 200 runs of a 10-fold CV 
#  rule: left (3) option with 1*MAD (see validation_sel_lambdas.R)
##
######################################################################
#
######################################################################
## PACKAGES AND SOURCE FILES
######################################################################
rm(list=ls())
library(pense)
library(parallel)
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(robustHD)
#setwd("CaseStudy/ValidationData")
Rcpp::sourceCpp('../tau_scale.cpp')
source("../aux_funcs.R")
source("../../simulation-study/ggplot-utils.R")
library("robust")
library("ROCR")
######################################################################
## DATA
######################################################################

#setwd("CaseStudy/ValidationData")
load("trainBiT1_testBiT2_a60.Rdata")

#TRAINING
x<-train 
rownames(x)<-x$sample
y<- x$`Max %DS LAD_clean`
x<- x%>% select(-c(sample,BIOBANK_ID,`Max %DS LAD_clean`))
x<-as.matrix(x)

##################################################################
## NO REGULARIZATION 
##################################################################
library("robust")
mm.all<-lmRob(y~x,efficiency = 0.85)
fit_vals.MM<-cbind(1, x) %*% coefficients(mm.all)
residuals.MM <- y-fit_vals.MM
resid_scale <- scale_tau(residuals.MM)
outs <- abs(residuals.MM) > 1.5 *resid_scale
aux_names <- rownames(outs)
aux_names[which(!outs)] <- rep("", length(which(!outs)))
aux_df <- data.frame(fit_vals = fit_vals.MM, resids = residuals.MM, outs = outs,
                     names = aux_names)

outliers.pensem <- aux_df %>% ggplot(aes(x = fit_vals, y = resids)) + geom_point() +
  geom_text(aes(label = names), hjust = 1.2) +
  geom_hline(yintercept = 0 ) +
  geom_hline(yintercept = 2 * resid_scale, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = -2 * resid_scale, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 1.5 * resid_scale, col = 'orange', linetype = 'dashed') +
  geom_hline(yintercept = -1.5 * resid_scale, col = 'orange', linetype = 'dashed') +
  scale_x_continuous(name = 'Fitted values') + scale_y_continuous(name = 'Residuals') #+

multi_save("../figures/outliers_MM_unweight_MRM.pdf", plot = outliers.pensem, width = 9.5, height = 3.5)

#TEST
x.new<-test
rownames(x.new)<-x.new$BIOBANK_ID
y.new<- x.new$`Max %DS LAD_clean`
x.new<- x.new%>% select(-c(sample,BIOBANK_ID,`Max %DS LAD_clean`))
x.new<-as.matrix(x.new)
pred_vals.MM <- cbind(1, x.new) %*% coefficients(mm.all)
residuals.test.MM <- y.new-pred_vals.MM
resid_scale <- scale_tau(residuals.test.MM)
outs <- abs(residuals.test.MM) > 1.5 *resid_scale
aux_names <- rownames(outs)
aux_names[which(!outs)] <- rep("", length(which(!outs)))
aux_df <- data.frame(fit_vals = pred_vals.MM, resids = residuals.test.MM, outs = outs,
                     names = aux_names)
outliers.pensem <- aux_df %>% ggplot(aes(x = fit_vals, y = resids)) + geom_point() +
  geom_text(aes(label = names), hjust = 1.2) +
  geom_hline(yintercept = 0 ) +
  geom_hline(yintercept = 2 * resid_scale, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = -2 * resid_scale, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 1.5 * resid_scale, col = 'orange', linetype = 'dashed') +
  geom_hline(yintercept = -1.5 * resid_scale, col = 'orange', linetype = 'dashed') +
  scale_x_continuous(name = 'Predicted values') + scale_y_continuous(name = 'Residuals') #+

multi_save("../figures/test_MM_unweight_MRM.pdf", plot = outliers.pensem, width = 9.5, height = 3.5)

test.clean<-test %>% filter(outs==0)
pred.clean<-pred_vals.MM[outs==0]
true_class<-1*(test.clean$`Max %DS LAD_clean`>30)
pred<-prediction(pred.clean,true_class)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

length(true_class)
sum(true_class)

# [[1]]
# [1] 0.8530466 (9 cases, 40 cases)
# Drops to 0.729798 if cutoff is 25 (18 cases, 23 controls). 

#permutation test
set.seed(123)
auc.random<-c()
for(i in 1:100){
  pred.clean.p<-sample(pred.clean)
  pred<-prediction(pred.clean.p,true_class)
  auc.perf.random = performance(pred, measure = "auc")
  auc.random<-c(auc.random,auc.perf.random@y.values)
}
# Permuted predictions give an AUC of .49
mean(as.vector(unlist(auc.random)))
#-----------------------------------------------------------------
# results from weighted PENSE(M)
# [[1]]
# [1] 0.8125 (10 cases, 33 cases)
# Drops to [1] 0.6366048 if cutoff is 20 (29 cases, 14 controls). 
# Drops to [1] 0.7159091 if cutoff is 25 (20 cases, 23 controls). 
# Permutation predictions can achieve an AUC of .57
#
######################################################################