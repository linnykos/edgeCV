rm(list=ls())
load("../results/vanishing_rank.RData")

# for each sparsity level, compute the amount overselected, exactly selected, and underselected for each n
mat_list <- lapply(1:4, function(i){
  tmp_res <- res[((i-1)*10+1):(i*10)]
  
  # loop over all n
  vec <- sapply(1:10, function(n){
    
    # clean results
    b <- sapply(tmp_res[[n]], function(x){all(is.na(x))})
    tmp_res2 <- tmp_res[[n]][!b]
    trials_mod <- length(tmp_res2)
    
    # loop over all the trials
    tmp_vec <- rep(0, 8)
    for(x in 1:trials_mod){
      # first do ecv
      idx <- which.min(tmp_res2[[x]]$err_vec)
      if(idx < 3) {tmp_vec[1] <- tmp_vec[1]+1
      } else if(idx > 3) {tmp_vec[3] <- tmp_vec[3]+1
      } else {tmp_vec[2] <- tmp_vec[2] + 1}
      
      # next do cvc
      idx <- which(tmp_res2[[x]]$p_vec <= 0.05)
      if(length(idx) == 0) {tmp_vec[7] <- tmp_vec[7]+1
      } else if(3 %in% idx) {tmp_vec[5] <- tmp_vec[5]+1
      } else if(all(idx < 3)){tmp_vec[4] <- tmp_vec[4]+1
      } else if(all(idx > 3)){tmp_vec[6] <- tmp_vec[6]+1
      } else tmp_vec[8] <- tmp_vec[8]+1
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

# plot the bars (from top to bottom: none, underest, exact ext, overest)
main_vec <- c("1", "1/n^0.25", "1/n^0.5", "log(n)/n")
for(i in 1:4){
  png(paste0("../figures/vanishing_rank_",i,".png"), height = 1200, width = 2000, res = 300, units = "px")
  par(mfrow = c(1,2), mar = c(4,4,4,0.5))
  
  # reformat and plot ecv
  tmp <- as.table(mat_list[[i]][1:3,])
  tmp <- tmp[c(3,2,1), ]
  rownames(tmp) <- c("Over","Exact","Under")
  colnames(tmp) <- as.character(paste0("n=",paramMat[1:10,1]))
  graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3)], main = paste0("Coherent-setting\nNCV: Rho = ", main_vec[i]))
  
  # reformat and plot cvc
  tmp <- mat_list[[i]][4:8,]
  tmp[4,] <- tmp[4,]+tmp[5,]
  tmp <- tmp[1:4,]
  tmp <- apply(tmp,2,function(x){x/sum(x)*50})
  tmp <- as.table(tmp)
  tmp <- tmp[c(3,2,1,4), ]
  rownames(tmp) <- c("Over","Exact","Under","None")
  colnames(tmp) <- as.character(paste0("n=",paramMat[1:10,1]))
  graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3,4)], main = paste0("Coherent-setting\nNCV+CVC: Rho = ", main_vec[i]))
  graphics.off()
}
