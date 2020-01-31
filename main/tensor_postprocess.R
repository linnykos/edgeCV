rm(list=ls())
load("../results/sparsity_5.RData")
alpha <- 0.05

k_vec <- as.character(3:5)
spar_level <- 6

# for each sparsity level, compute the amount overselected, exactly selected, and underselected for each n
mat_list <- lapply(1:length(k_vec), function(i){
  tmp_res <- res[((i-1)*spar_level+1):(i*spar_level)]
  
  
  vec <- sapply(1:length(tmp_res), function(n){
    # clean results
    b <- sapply(tmp_res[[n]], function(x){all(is.na(x))})
    tmp_res2 <- tmp_res[[n]][!b]
    trials_mod <- length(tmp_res2)
    
    # loop over all the trials
    tmp_vec <- rep(0, 12)
    for(x in 1:trials_mod){
      # first do ecv
      idx <- which.min(tmp_res2[[x]]$ecv_err_vec)
      if(idx < k_vec[i]) {tmp_vec[1] <- tmp_vec[1]+1
      } else if(idx > k_vec[i]) {tmp_vec[3] <- tmp_vec[3]+1
      } else {tmp_vec[2] <- tmp_vec[2] + 1}
      
      # next do cvc
      idx <- which(tmp_res2[[x]]$ecv_p_vec >= alpha)
      if(length(idx) == 0) {
        tmp_vec[7] <- tmp_vec[7]+1
      } else if(k_vec[i] %in% idx) {tmp_vec[5] <- tmp_vec[5]+1
      } else if(all(idx < k_vec[i])){tmp_vec[4] <- tmp_vec[4]+1
      } else if(all(idx > k_vec[i])){tmp_vec[6] <- tmp_vec[6]+1
      } else tmp_vec[7] <- tmp_vec[7] + 1 #what to do when results are logically incoherent
      
      idx <- idx[which.min(idx)]
      if(length(idx) == 0) {
        tmp_vec[12] <- tmp_vec[12]+1
      } else if(idx ==  k_vec[i]) {tmp_vec[10] <- tmp_vec[10]+1
      } else if(all(idx <  k_vec[i])){tmp_vec[9] <- tmp_vec[9]+1
      } else if(all(idx >  k_vec[i])){tmp_vec[11] <- tmp_vec[11]+1
      }
    }
    
    tmp_vec
  })
  
  vec
})

col_vec <- c(rgb(165, 217, 151, maxColorValue = 255), #green
             rgb(117, 164, 242, maxColorValue = 255), #blue
             rgb(239, 133, 140, maxColorValue = 255), #red
             rgb(73, 73, 73, maxColorValue = 255) #gray
)