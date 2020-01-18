rm(list=ls())

set.seed(10)
b_mat_truth <- matrix(c(1/4, 1/2, 1/4,  1/2, 1/4, 1/4,  1/4, 1/4, 1/6), 3, 3)
cluster_idx_truth <- sample(c(1:3), 150, prob = c(0.3, 0.3, 0.4), replace = T)
dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
k <- 5
ncv_res <- network_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)
err_mat <- ncv_res$err_mat
trials <- 100
alpha <- 0.05
verbose = T

############

stopifnot(length(colnames(err_mat)) == ncol(err_mat))

n <- nrow(err_mat)
k <- ncol(err_mat)
err_mat2 <- .clean_err_mat(err_mat)

mu_vec <- colMeans(err_mat2)
sd_vec <- apply(err_mat2, 2, sd)
test_vec <- .extract_max(sqrt(n)*mu_vec/sd_vec, k)

boot_mat <- sapply(1:trials, function(b){
  if(verbose && b %% floor(trials/10) == 0) cat('*')
  .cvc_bootstrap_trial(err_mat2, mu_vec, sd_vec, k)
})