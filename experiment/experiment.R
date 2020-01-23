rm(list=ls())
set.seed(10)
b_mat_truth <- diag(3)
cluster_idx_truth <- rep(1:3, each = 50)
dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
k <- 5
err_mat_list <- edge_cv_sbm(dat, k_vec = c(1:k), nfold = 5, verbose = F)$err_mat_list
trials <- 100
alpha <- 0.05
verbose = T

for(i in 1:length(err_mat_list)){
  stopifnot(length(colnames(err_mat_list[[i]])) == ncol(err_mat_list[[i]]))
}

n_vec <- sapply(err_mat_list, function(err_mat){nrow(err_mat)})
k_vec <- sapply(err_mat_list, function(err_mat){ncol(err_mat)})
stopifnot(length(unique(k_vec)) == 1)
k <- unique(k_vec)
n <- sum(sapply(err_mat_list, nrow))

err_mat_list2 <- lapply(err_mat_list, function(err_mat){
  .clean_err_mat(err_mat)
})

mu_vec_list <- lapply(err_mat_list2, function(err_mat2){colMeans(err_mat2)})

err_mat_list2_recentered <- .recenter_list(err_mat_list2, mu_vec_list)

mu_vec <- colMeans(do.call(rbind, mu_vec_list))
sd_vec <- apply(do.call(rbind, err_mat_list2_recentered), 2, stats::sd)

test_vec <- .extract_max(sqrt(n)*mu_vec/sd_vec, k)

boot_mat <- sapply(1:trials, function(b){
  if(verbose && b %% floor(trials/10) == 0) cat('*')
  
  tmp_mat <- .cvc_bootstrap_trial_cv(err_mat_list2, mu_vec_list, sd_vec)
  .extract_max((1/sqrt(n)) * colSums(tmp_mat), k)
})

p_vec <- sapply(1:k, function(x){
  length(which(boot_mat[x,] <= test_vec[x]))/ncol(boot_mat)
})