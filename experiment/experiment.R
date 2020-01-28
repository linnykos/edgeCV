set.seed(10)
beta <- 0.25
b_mat_truth <- (1-1.5*beta)*diag(3) + beta*rep(1,3) %*% t(rep(1,3))
cluster_idx_truth <- rep(1:3, each = 10)
dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
k <- 5
ecv_res <- edge_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)
err_mat <- ecv_res$err_mat

B = 200
alpha = 0.05

###########

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

sgmb_p_val
