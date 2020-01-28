rm(list=ls())
load("../results/tensor.RData")

n_level <- 3
spar_level <- 3

mat_list <- lapply(1:spar_level, function(i){
  tmp_res <- res[((i-1)*n_level+1):(i*n_level)]
  
  # loop over all n
  vec <- sapply(1:length(tmp_res), function(n){
    
    # clean results
    b <- sapply(tmp_res[[n]], function(x){all(is.na(x))})
    tmp_res2 <- tmp_res[[n]][!b]
    trials_mod <- length(tmp_res2)
    
    # loop over all the trials
    tmp_vec <- rep(0, 3)
    for(x in 1:trials_mod){
      # first do ecv
      idx <- which.min(tmp_res2[[x]]$err_vec)
      if(idx < 3) {tmp_vec[1] <- tmp_vec[1]+1
      } else if(idx > 3) {tmp_vec[3] <- tmp_vec[3]+1
      } else {tmp_vec[2] <- tmp_vec[2] + 1}
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

main_vec <- c("1", "1/n^0.25", "1/n^0.5")
n_vec <- as.character(c(paramMat[1:3]))
png("../figures/tensor.png", height = 1200, width = 3000, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
for(i in 1:spar_level){
  
  # reformat and plot ecv
  tmp <- as.table(mat_list[[i]][1:3,])
  tmp <- apply(tmp,2,function(x){x/sum(x)*50})
  tmp <- tmp[c(3,2,1), ]
  rownames(tmp) <- c("Over","Exact","Under")
  colnames(tmp) <- as.character(paste0("n=",n_vec[1:ncol(mat_list[[i]])]))
  graphics::barplot(as.table(tmp), horiz = TRUE, col = col_vec[c(2,1,3)], 
                    main = paste0("Tensor SBM-setting\nECV: Rho = ", main_vec[i]))
}
graphics.off()
