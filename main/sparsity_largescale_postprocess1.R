rm(list=ls())
prefix <- "../results/sparsity_largescale_"
counter <- 1:10

res_all <- vector("list", 1)
res_all[[1]] <- vector("list", 1)
current_len <- 0

for(counter_idx in counter){
  load(paste0(prefix, counter_idx, ".RData"))
  len <- length(res_list)
  
  res_all[[1]][(current_len+1):(current_len+len)] <- res_list
  
  current_len <- length(res_all[[1]])
}

res <- res_all
rm(list=c("res_all", "res_list"))
save.image("sparsity_largescale_all.RData")