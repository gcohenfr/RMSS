#---------------------------------------
# GENERATE FIGURES FROM THE SIMULATIONS
#---------------------------------------

# Clear Memory
rm(list=ls())

# Required Library
library(ggplot2)
library(tidyverse)
library(Cairo)
library(ggh4x) 

# Required Source Files


#---------------------------------------

# __________
# MSPE Plot
# __________

# Plot data (MSPE) - No Contamination
method.names <- c("PENSE", "RBSS", "RMSS")
plot.data <- matrix(nrow = 0, ncol = 5)
load("results/results_n=50_p=500_cont=Cohen-Freue_snr=1_rho=0.8_tau=0.Rdata")
for(iter in 1:50){
  
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.1", "No Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.2", "No Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.4", "No Contamination"))
}
# Plot data (MSPE) - Moderate Contamination
load("results/results_n=50_p=500_cont=Cohen-Freue_snr=1_rho=0.8_tau=0.15.Rdata")
for(iter in 1:50){
  
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.1", "Moderate Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.2", "Moderate Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.4", "Moderate Contamination"))
}
# Plot data (MSPE) - High Contamination
load("results/results_n=50_p=500_cont=Cohen-Freue_snr=1_rho=0.8_tau=0.3.Rdata")
for(iter in 1:50){
  
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.1", "High Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.2", "High Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "MSPE", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "MSPE"], "0.4", "High Contamination"))
}

# Plot metadata
plot.data <- data.frame(plot.data)
rownames(plot.data) <- NULL
colnames(plot.data) <- c("Method", "Measure", "MSPE", "Sparsity_Level", "Contamination_Category")
plot.data$Method <- factor(plot.data$Method,
                           levels = c("PENSE", "RBSS", "RMSS"), 
                           ordered = TRUE)
plot.data$Measure <- factor(plot.data$Measure, 
                                   levels=c("MSPE"))
plot.data$Sparsity_Level <- factor(plot.data$Sparsity_Level, 
                                   levels=c("0.1", "0.2", "0.4"))
plot.data$Contamination_Category <- factor(plot.data$Contamination_Category, 
                                           levels=c("No Contamination", "Moderate Contamination", "High Contamination"))
plot.data$MSPE <- as.numeric(as.character(plot.data$MSPE))

# Remove extreme MSPEs
plot.data <- plot.data[!(plot.data$Contamination_Category == "No Contamination" & plot.data$MSPE > 2),]
plot.data <- plot.data[!(plot.data$Contamination_Category == "Moderate Contamination" & plot.data$MSPE > 2.5),]
plot.data <- plot.data[!(plot.data$Contamination_Category == "High Contamination" & plot.data$MSPE > 7.5),]

# Plot 
plot.data %>% 
  ggplot(aes(x=Sparsity_Level, y=MSPE, aes=Contamination_Category, fill=Method)) +
  stat_boxplot(geom="errorbar", position = position_dodge(width=0.5), width = 0.35) +
  geom_boxplot(position = position_dodge(width=0.5), width=0.35) + 
  facet_grid2(Measure ~ Contamination_Category, scales = "free_y", independent = "y") +
  labs(fill = "Method", x="") + 
  scale_fill_manual(values=c("#696969", "#979797", "#C5C5C5")) + 
  theme_bw() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  theme(text = element_text(size = 12)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sparsity Level") + 
  ylab("")

# Saving plot
ggsave("Simulation_Plot_MSPE.png", device = "png", dpi=600)

# ___________
# RC/PR Plot
# ___________

# Plot data (RC/PR) - No Contamination
method.names <- c("PENSE", "RBSS", "RMSS")
plot.data <- matrix(nrow = 0, ncol = 5)
load("results/results_n=50_p=500_cont=Cohen-Freue_snr=1_rho=0.8_tau=0.Rdata")
for(iter in 1:50){
  
  # RC data
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.1", "No Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.2", "No Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.4", "No Contamination"))
  
  # PR data
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.1", "No Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.2", "No Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.4", "No Contamination"))
}
# Plot data (RC/PR) - Moderate Contamination
load("results/results_n=50_p=500_cont=Cohen-Freue_snr=1_rho=0.8_tau=0.15.Rdata")
for(iter in 1:50){
  
  # RC data
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.1", "Moderate Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.2", "Moderate Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.4", "Moderate Contamination"))
  
  # PR data
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.1", "Moderate Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.2", "Moderate Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.4", "Moderate Contamination"))
}
# Plot data (RC/PR) - High Contamination
load("results/results_n=50_p=500_cont=Cohen-Freue_snr=1_rho=0.8_tau=0.3.Rdata")
for(iter in 1:50){
  
  # RC data
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.1", "High Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.2", "High Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "RC", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "RC"], "0.4", "High Contamination"))
  
  # PR data
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[1]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.1", "High Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[2]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.2", "High Contamination"))
  plot.data <- rbind(plot.data, 
                     cbind(method.names, "PR", output[[3]][,, iter][c("PENSE", "RBSS", "RMSS"), "PR"], "0.4", "High Contamination"))
  
}

# Plot metadata
plot.data <- data.frame(plot.data)
rownames(plot.data) <- NULL
colnames(plot.data) <- c("Method", "Measure", "Data", "Sparsity_Level", "Contamination_Category")
plot.data$Method <- factor(plot.data$Method,
                           levels = c("PENSE", "RBSS", "RMSS"), 
                           ordered = TRUE)
plot.data$Measure <- factor(plot.data$Measure,
                            levels = c("RC", "PR"),
                            ordered = TRUE)
plot.data$Sparsity_Level <- factor(plot.data$Sparsity_Level, 
                                   levels=c("0.1", "0.2", "0.4"))
plot.data$Contamination_Category <- factor(plot.data$Contamination_Category, 
                                           levels=c("No Contamination", "Moderate Contamination", "High Contamination"))
plot.data$Data <- as.numeric(as.character(plot.data$Data))

# Plot 
plot.data %>% 
  ggplot(aes(x=Sparsity_Level, y=Data, aes=Contamination_Category, fill=Method)) +
  stat_boxplot(geom="errorbar", position = position_dodge(width=0.5), width = 0.35) +
  geom_boxplot(position = position_dodge(width=0.5), width=0.35) + 
  facet_grid2(Measure ~ Contamination_Category, scales = "free_y", independent = "y") +
  labs(fill = "Method", x="") + 
  scale_fill_manual(values=c("#696969", "#979797", "#C5C5C5")) + 
  theme_bw() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  theme(text = element_text(size = 12)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Sparsity Level") + 
  ylab("")

# Saving plot
ggsave("Simulation_Plot_RCPR.png", device = "png", dpi=600)

