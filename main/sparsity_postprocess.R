rm(list=ls())
load("../results/sparsity_4.RData")
alpha <- 0.05

k_level <- 3
spar_level <- 11

# for each sparsity level, compute the amount overselected, exactly selected, and underselected for each n
mat_list <- lapply(1:k_level, function(i){
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
      idx <- which.min(tmp_res2[[x]]$err_vec)
      if(idx < 3) {tmp_vec[1] <- tmp_vec[1]+1
      } else if(idx > 3) {tmp_vec[3] <- tmp_vec[3]+1
      } else {tmp_vec[2] <- tmp_vec[2] + 1}
      
      # next do cvc
      idx <- which(tmp_res2[[x]]$p_vec >= alpha)
      if(length(idx) == 0) {
        tmp_vec[7] <- tmp_vec[7]+1
      } else if(3 %in% idx) {tmp_vec[5] <- tmp_vec[5]+1
      } else if(all(idx < 3)){tmp_vec[4] <- tmp_vec[4]+1
      } else if(all(idx > 3)){tmp_vec[6] <- tmp_vec[6]+1
      } else tmp_vec[7] <- tmp_vec[7] + 1 #what to do when results are logically incoherent
      
      idx <- idx[which.min(idx)]
      if(length(idx) == 0) {
        tmp_vec[12] <- tmp_vec[12]+1
      } else if(idx == 3) {tmp_vec[10] <- tmp_vec[10]+1
      } else if(all(idx < 3)){tmp_vec[9] <- tmp_vec[9]+1
      } else if(all(idx > 3)){tmp_vec[11] <- tmp_vec[11]+1
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

# plot the bars (from top to bottom: none, underest, exact ext, overest)
main_vec <- c("1", paste0("1/n^", paramMat[2:11,2]))
k_vec <- as.character(3:5)
for(k in 1:k_level){
  png(paste0("../figures/sparsity_4_", k, ".png"), height = 1500, width = 4500, res = 300, units = "px")
  par(mfrow = c(1,3), mar = c(4,4,4,1))
  
  # reformat and plot ecv
  tmp <- mat_list[[k]][1:3,]
  tmp <- apply(tmp,2,function(x){x/sum(x)*trials})
  tmp <- tmp[c(3,2,1), ]
  rownames(tmp) <- c("Over","Exact","Under")
  colnames(tmp) <- main_vec
  graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3)], 
                    main = paste0("SBM-setting, ECV\n", k_vec[k], " clusters"),
                    cex.names = 0.7)
  
  # reformat and plot cvc
  tmp <- mat_list[[k]][4:8,]
  tmp <- apply(tmp,2,function(x){x/sum(x)*trials})
  tmp <- tmp[c(3,2,1,4), ]
  rownames(tmp) <- c("Over","Exact","Under","None")
  colnames(tmp) <- main_vec
  graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3,4)], 
                    main = paste0("SBM-setting, ECV+CVC\n", k_vec[k], " clusters"),
                    cex.names = 0.7)
  
  # reformat and plot cvc
  tmp <- mat_list[[k]][9:12,]
  tmp <- apply(tmp,2,function(x){x/sum(x)*trials})
  tmp <- tmp[c(3,2,1,4), ]
  rownames(tmp) <- c("Over","Exact","Under","None")
  colnames(tmp) <- main_vec
  graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3,4)], 
                    main = paste0("SBM-setting, ECV+CVC\n", k_vec[k], " clusters"),
                    cex.names = 0.7)
  graphics.off()
}

########

success_idx <- c(2,5,10)
png(paste0("../figures/sparsity_byn.png"), height = 1500, width = 4500, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,4,1))
for(k in 1:k_level){
  plot(NA, xlim = range(paramMat[1:11,2]), ylim = c(0, 1), xlab = "Sparsity level: (1/n^x)", ylab = "Percentage success",
       main = paste0("SBM with ", k_vec[k], " clusters"))
  for(i in 1:3){
    lines(x = paramMat[1:11,2], y = mat_list[[k]][success_idx[i],]/trials, col = i, lwd = 2)
    points(x = paramMat[1:11,2], y = mat_list[[k]][success_idx[i],]/trials, col = i, pch = 16, cex = 2)
  }
}
