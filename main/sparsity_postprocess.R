rm(list=ls())
load("../results/sparsity_3.RData")
alpha <- 0.05

# for each sparsity level, compute the amount overselected, exactly selected, and underselected for each n
mat_list <- lapply(1:nrow(paramMat), function(i){
  tmp_res <- res[[i]]
  
  # clean results
  b <- sapply(tmp_res, function(x){all(is.na(x))})
  tmp_res2 <- tmp_res[!b]
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
    idx <- which(tmp_res2[[x]]$p_vec >= alpha)
    if(length(idx) == 0) {
      tmp_vec[4] <- tmp_vec[4]+1
    } else if(3 %in% idx) {tmp_vec[5] <- tmp_vec[5]+1
    } else if(all(idx < 3)){tmp_vec[4] <- tmp_vec[4]+1
    } else if(all(idx > 3)){tmp_vec[6] <- tmp_vec[6]+1
    } else tmp_vec[4] <- tmp_vec[4]+1 #what to do when results are logically incoherent
    # idx <- idx[which.min(idx)]
    # if(length(idx) == 0) {
    #   tmp_vec[7] <- tmp_vec[7]+1
    # } else if(idx == 3) {tmp_vec[5] <- tmp_vec[5]+1
    # } else if(all(idx < 3)){tmp_vec[4] <- tmp_vec[4]+1
    # } else if(all(idx > 3)){tmp_vec[6] <- tmp_vec[6]+1
    # }
  }
  
  tmp_vec
})

col_vec <- c(rgb(165, 217, 151, maxColorValue = 255), #green
             rgb(117, 164, 242, maxColorValue = 255), #blue
             rgb(239, 133, 140, maxColorValue = 255), #red
             rgb(73, 73, 73, maxColorValue = 255) #gray
             )

res_mat <- do.call(rbind, mat_list)

# plot the bars (from top to bottom: none, underest, exact ext, overest)
main_vec <- c("1", paste0("1/n^", paramMat[2:11,2]), "log(n)/n")
png(paste0("../figures/sparsity_3.png"), height = 1500, width = 3000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,4,1))

# reformat and plot ecv
tmp <- t(res_mat[,1:3])
tmp <- apply(tmp,2,function(x){x/sum(x)*trials})
tmp <- tmp[c(3,2,1), ]
rownames(tmp) <- c("Over","Exact","Under")
colnames(tmp) <- main_vec
graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3)], 
                  main = paste0("SBM-setting, ECV"),
                  cex.names = 0.7)

# reformat and plot cvc
tmp <- t(res_mat[,4:8])
tmp <- apply(tmp,2,function(x){x/sum(x)*trials})
tmp <- tmp[c(3,2,1,4), ]
rownames(tmp) <- c("Over","Exact","Under","None")
colnames(tmp) <- main_vec
graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3,4)], 
                  main = paste0("SBM-setting, ECV+CVC"),
                  cex.names = 0.7)
graphics.off()

