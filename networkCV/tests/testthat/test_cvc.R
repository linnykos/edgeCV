context("Test cvc")

## cvc is correct

test_that("cvc_sbm works", {
  set.seed(10)
  b_mat_truth <- matrix(c(1/4, 1/2, 1/4,  1/2, 1/4, 1/4,  1/4, 1/4, 1/6), 3, 3)
  cluster_idx_truth <- c(rep(1, 10), rep(2, 10), rep(3, 10))
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  err_mat_list <- edge_cv_sbm(dat, k_vec = c(1:k), nfold = 5, verbose = F)$err_mat_list
  trials <- 10
  alpha <- 0.05
  
  res <- cvc(do.call(rbind, err_mat_list))
  
  expect_true(length(res) == 5)
  expect_true(all(is.numeric(res)))
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
})

test_that("cvc_sbm works with more samples", {
  set.seed(10)
  b_mat_truth <- matrix(c(1/4, 1/2, 1/4,  1/2, 1/4, 1/4,  1/4, 1/4, 1/6), 3, 3)
  n_each <- 200
  cluster_idx_truth <- c(rep(1, n_each), rep(2, n_each), rep(3, n_each))
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  err_mat_list <- edge_cv_sbm(dat, k_vec = c(1:k), nfold = 5, verbose = F)$err_mat_list
  alpha <- 0.05
  
  res <- cvc(do.call(rbind, err_mat_list))
  
  expect_true(length(res) == 5)
  expect_true(all(is.numeric(res)))
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
})

test_that("cvc_sbm works can havef non-trivial p-values", {
  set.seed(10)
  b_mat_truth <- matrix(c(1/4, 1/2, 1/4,  1/2, 1/4, 1/4,  1/4, 1/4, 1/6), 3, 3)
  n_each <- 5
  cluster_idx_truth <- c(rep(1, n_each), rep(2, n_each), rep(3, n_each))
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  k <- 5
  err_mat_list <- edge_cv_sbm(dat, k_vec = c(1:k), nfold = 5, verbose = F)$err_mat_list
  alpha <- 0.05
  
  res <- cvc(do.call(rbind, err_mat_list))
  
  expect_true(all(res > 0))
  expect_true(all(res < 1))
})