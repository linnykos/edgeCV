context("Test cvc")

## .clean_err_mat is correct

test_that(".clean_err_mat works", {
  set.seed(10)
  b_mat_truth <- 0.5*diag(3)
  b_mat_truth <- b_mat_truth + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  ncv_output <- network_cv_sbm_sample_split(dat, k_vec = c(1:5), test_prop = 0.1)
  
  res <- .clean_err_mat(ncv_output$err_mat)
  
  expect_true(nrow(res) == nrow(ncv_output$err_mat))
  k <- ncol(ncv_output$err_mat)
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
  ncv_output <- network_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)
  
  err_mat2 <- .clean_err_mat(ncv_output$err_mat)
  
  mu_vec <- colMeans(err_mat2)
  sd_vec <- apply(err_mat2, 2, sd)
  
  res <- .cvc_bootstrap_trial(err_mat2, mu_vec, sd_vec, k)
  
  expect_true(length(res) == k)
})

#################

## cvc_sbm_sample_split is correct

test_that("cvc_sbm_sample_split works", {
  set.seed(10)
  b_mat_truth <- diag(3)
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  err_mat <- network_cv_sbm_sample_split(dat, k_vec = c(1:k), test_prop = 0.1)$err_mat
  trials <- 200
  alpha <- 0.5
  
  res <- cvc_sbm_sample_split(err_mat, trials, alpha)
  
  expect_true(is.list(res))
  expect_true(all(is.character(res$model_selected)))
  expect_true(length(res$p_vec) == 5)
})