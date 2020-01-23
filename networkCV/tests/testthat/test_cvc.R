context("Test cvc")

## .clean_err_mat is correct

test_that(".clean_err_mat works", {
  set.seed(10)
  b_mat_truth <- 0.5*diag(3)
  b_mat_truth <- b_mat_truth + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  ecv_output <- edge_cv_sbm_sample_split(dat, k_vec = c(1:5), test_prop = 0.1)
  
  res <- .clean_err_mat(ecv_output$err_mat)
  
  expect_true(nrow(res) == nrow(ecv_output$err_mat))
  k <- ncol(ecv_output$err_mat)
  expect_true(ncol(res) == k*(k-1))
})


################

## .extract_max is correct

test_that(".extract_max works", {
  vec <- 1:90
  res <- .extract_max(vec, 10)
  
  expect_true(length(res) == 10)
  expect_true(sum(is.na(res)) == 0)
})


##############3

## .cvc_bootstrap_trial is correct

test_that(".cvc_bootstrap_trial works", {
  set.seed(10)
  b_mat_truth <- 0.5*diag(3)
  b_mat_truth <- b_mat_truth + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  ecv_output <- edge_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)
  
  err_mat2 <- .clean_err_mat(ecv_output$err_mat)
  
  mu_vec <- colMeans(err_mat2)
  sd_vec <- apply(err_mat2, 2, sd)
  
  res <- .cvc_bootstrap_trial(err_mat2, mu_vec, sd_vec, k)
  
  expect_true(all(dim(res) == dim(err_mat2)))
})

test_that(".cvc_bootstrap_trial makes the correct calculation", {
  set.seed(10)
  b_mat_truth <- 0.5*diag(3)
  b_mat_truth <- b_mat_truth + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  ecv_output <- edge_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)
  
  err_mat2 <- .clean_err_mat(ecv_output$err_mat)
  
  mu_vec <- colMeans(err_mat2)
  sd_vec <- apply(err_mat2, 2, sd)
  
  set.seed(10)
  res <- .cvc_bootstrap_trial(err_mat2, mu_vec, sd_vec, k)
  
  set.seed(10)
  g_vec <- stats::rnorm(nrow(err_mat2))
  tmp_mat <- matrix(NA, nrow = nrow(err_mat2), ncol = ncol(err_mat2))
  for(i in 1:nrow(err_mat2)){
    for(j in 1:ncol(err_mat2)){
      tmp_mat[i,j] <- (err_mat2[i,j] - mu_vec[j])*g_vec[i]/sd_vec[j]
    }
  }
  
  expect_true(sum(abs(res-  tmp_mat)) <= 1e-6)
})

#################

## cvc_sbm_sample_split is correct

test_that("cvc_sbm_sample_split works", {
  set.seed(10)
  b_mat_truth <- diag(3)
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  err_mat <- edge_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)$err_mat
  trials <- 100
  alpha <- 0.05
  
  res <- cvc_sbm_sample_split(err_mat, trials, alpha)
  
  expect_true(is.list(res))
  expect_true(all(is.character(res$model_selected)))
  expect_true(length(res$p_vec) == 5)
})

test_that("cvc_sbm_sample_split works on a slightly more interesting example", {
  set.seed(10)
  beta <- 0.25
  b_mat_truth <- (1-1.5*beta)*diag(3) + beta*rep(1,3) %*% t(rep(1,3))
  cluster_idx_truth <- c(rep(1, 45), rep(2, 45), rep(3, 60))
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  err_mat <- edge_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)$err_mat
  trials <- 100
  alpha <- 0.05
  
  res <- cvc_sbm_sample_split(err_mat, trials, alpha)
  
  expect_true(is.list(res))
  expect_true(all(is.character(res$model_selected)))
  expect_true(length(res$p_vec) == 5)
})

test_that("cvc_sbm_sample_split works when the difference between ranks is small", {
  set.seed(10)
  b_mat_truth <- matrix(c(0.8, 0.7, 0.1,  0.7, 0.8, 0.1, 0.1, 0.1, 0.8), 3, 3)
  cluster_idx_truth <- c(rep(1, 45), rep(2, 45), rep(3, 60))
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  ecv_res <- edge_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)
  err_mat <- ecv_res$err_mat
  trials <- 200
  alpha <- 0.05
  
  res <- cvc_sbm_sample_split(err_mat, trials, alpha)
  
  expect_true(is.list(res))
  expect_true(all(is.character(res$model_selected)))
  expect_true(length(res$p_vec) == 5)
})

test_that("cvc_sbm_sample_split works with negative eigenvalues", {
  set.seed(10)
  b_mat_truth <- matrix(c(1/4, 1/2, 1/4,  1/2, 1/4, 1/4,  1/4, 1/4, 1/6), 3, 3)
  cluster_idx_truth <- c(rep(1, 45), rep(2, 45), rep(3, 60))
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  err_mat <- edge_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)$err_mat
  trials <- 100
  alpha <- 0.05
  
  res <- cvc_sbm_sample_split(err_mat, trials, alpha)
  
  expect_true(is.list(res))
  expect_true(all(is.character(res$model_selected)))
  expect_true(length(res$p_vec) == 5)
})
