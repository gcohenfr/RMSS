#--------------------------------
# VARIABLE SELECTION - BBS DATA
#--------------------------------

# Clear Memory
rm(list=ls())

# Required Library


# Required Source Files


#--------------------------------

# RMSS - Genes most selected
rmss_selections <- numeric(500)
for(iter in 1:50){
  
  rmss_selections[as.numeric(rownames(selections$RMSS[[iter]]))] <- 
    rmss_selections[as.numeric(rownames(selections$RMSS[[iter]]))] + 1
}
rmss_top <- order(rmss_selections, decreasing = TRUE)
rmss_top[1:50]
rmss_selections[rmss_top[1:50]]

# PENSE - Genes most selected
pense_selections <- numeric(500)
for(iter in 1:50){
  
  pense_selections[selections$PENSE[[iter]]] <- 
    pense_selections[selections$PENSE[[iter]]] + 1
}
pense_top <- order(pense_selections, decreasing = TRUE)
pense_top[1:50]
pense_selections[pense_top[1:50]]

# Sparse LTS - Genes most selected
sparselts_selections <- numeric(500)
for(iter in 1:50){
  
  sparselts_selections[selections$SparseLTS[[iter]]] <- 
    sparselts_selections[selections$SparseLTS[[iter]]] + 1
}
sparselts_top <- order(sparselts_selections, decreasing = TRUE)
sparselts_top[1:50]
sparselts_selections[sparselts_top[1:50]]

# Comparing top genes
pense_selections[pense_top[1:50]]
pense_selections[rmss_top[1:50]]

sparselts_selections[sparselts_top[1:50]]
sparselts_selections[rmss_top[1:50]]

rmss_selections[rmss_top[1:50]]
rmss_selections[pense_top[1:50]]
rmss_selections[sparselts_top[1:50]]


