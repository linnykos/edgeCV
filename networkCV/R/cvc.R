cvc <- function(err_mat, B = 200, alpha = 0.05) {
  clean_res <- .clean_err_mat(err_mat)
  
  err_mat_c <- clean_res$err_mat_c
  err_mat_c_ind <- clean_res$ind
  n <- nrow(err_mat_c)
  M <- ncol(err_mat_c)
  
  err_mean <- apply(err_mat_c, 2, mean)
  err_mat_center <- err_mat_c - matrix(err_mean, nrow = n, ncol = M, byrow = T)
  
  gauss_mat <- matrix(stats::rnorm(n*B), ncol = B)
  
  sgmb_p_val <- sapply(1:M, function(m){
    err_diff_center <- err_mat_center[,m] - err_mat_center[,-m, drop=F]
    err_mean_diff <- err_mean[m] - err_mean[-m]
    
    sd_vec <- apply(err_diff_center, 2, sd)
    err_mean_diff_scale <- sqrt(n) * err_mean_diff / sd_vec
    
    err_diff_center_scale <- err_diff_center / 
      matrix(sd_vec, nrow = n, ncol = ncol(err_diff_center), byrow = T)
    
    # compute the bootstrap statistic
    test_stat_vec <- sapply(1:B, function(ib){
      sqrt(n) * max(apply(err_diff_center_scale*gauss_mat[,ib], 2, mean))
    })
    
    mean(test_stat_vec > max(err_mean_diff_scale))
  })
  
  res <- rep(0, ncol(err_mat))
  for (i in 1:length(err_mat_c_ind)) {
    res[err_mat_c_ind[[i]]] <- sgmb_p_val[i]
  }
  
  res
}


##############

.clean_err_mat <- function(err_mat) {
  n2 <- nrow(err_mat)
  M <- ncol(err_mat)
  
  # center the matrix
  err_mean <- apply(err_mat, 2, mean)
  err_mat_center <- err_mat - matrix(err_mean, nrow = n2, ncol = M, byrow = T)
  
  res <- list()
  list_len <- 0
  remaining_index <- 1:M
  
  while (length(remaining_index) > 0) {
    m <- remaining_index[1]
    err_diff_center <- err_mat_center[,m] - err_mat_center
    err_mean_diff <- err_mean[m] - err_mean
    sd_vec <- apply(err_diff_center, 2, sd)
    J_1 <- which(sd_vec==0 & err_mean_diff > 0)
    J_2 <- which(sd_vec==0 & err_mean_diff == 0)
    J_3 <- which(sd_vec==0 & err_mean_diff < 0)
    
    if (length(J_1) == 0) {
      list_len <- list_len+1
      res[[list_len]] <- J_2
    }
    
    remaining_index <- setdiff(remaining_index, union(J_2, J_3))
  }
  
  list(err_mat_c = err_mat[,sapply(res,function(x){x[1]})], ind=res)
}