context("Test network cv")

## .remove_entries is correct

test_that(".remove_entries works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx <- rep(1:3, each = 5)
  dat <- generate_sbm(b_mat, cluster_idx)
  
  res <- .remove_entries(dat, 0.5)
  
  expect_true(all(dim(dat) == dim(res)))
})

test_that(".remove_entries removes roughly the correct amount of entries", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx <- rep(1:3, each = 5)
  dat <- generate_sbm(b_mat, cluster_idx)
  n <- nrow(dat)
  
  res <- .remove_entries(dat, 0.4)
  
  num_NA <- length(which(is.na(res)))
  
  expect_true(abs(num_NA - 2*0.4*choose(n,2)) <= 1)
})

################

## .impute_matrix is correct

test_that(".impute_matrix works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat, cluster_idx)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  
  res <- .impute_matrix(dat_NA, ncol(b_mat), test_prop)
  
  expect_true(all(dim(dat) == dim(res)))
})

test_that(".impute_matrix works reasonably well", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat, cluster_idx)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  res <- .impute_matrix(dat_NA, ncol(b_mat), test_prop)
  
  b_mat2 <- matrix(0.7, 3, 3)
  b_mat2 <- b_mat2 - 0.5*diag(3)
  dat2 <- generate_sbm(b_mat2, cluster_idx)
  dat_NA2 <- .remove_entries(dat2, test_prop)
  res2 <- .impute_matrix(dat_NA2, ncol(b_mat2), test_prop)
  
  cluster_mat <- .convert_idx_to_mat(cluster_idx)
  pop_mat <- cluster_mat %*% b_mat %*% t(cluster_mat)
  
  expect_true(sum(abs(pop_mat - res)) <= sum(abs(pop_mat - res2)))
})

################

## .spectral_clustering is correct

test_that(".spectral_clustering works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat, cluster_idx)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  dat_impute <- .impute_matrix(dat_NA, ncol(b_mat), test_prop)
  
  res <- .spectral_clustering(dat_impute, ncol(b_mat))
  
  expect_true(length(res) == ncol(dat))
  expect_true(all(res %% 1 == 0))
  expect_true(all(res >= 1))
  expect_true(all(res <= ncol(b_mat)))
})

test_that(".spectral_clustering can give correct answer", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat, cluster_idx)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  dat_impute <- .impute_matrix(dat_NA, ncol(b_mat), test_prop)
  
  res <- .spectral_clustering(dat_impute, ncol(b_mat))
  
  expect_true(length(unique(res[1:50])) == 1)
  expect_true(length(unique(res[51:100])) == 1)
  expect_true(length(unique(res[101:150])) == 1)
  expect_true(length(unique(res)) == 3)
})

#####################

## .form_prediction_sbm works

test_that(".form_prediction_sbm works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat, cluster_idx_truth)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  dat_impute <- .impute_matrix(dat_NA, ncol(b_mat), test_prop)
  cluster_idx <- .spectral_clustering(dat_impute, ncol(b_mat))
  
  res <- .form_prediction_sbm(dat_impute, cluster_idx)
  
  expect_true(all(dim(res) == dim(b_mat)))
})

test_that(".form_prediction_sbm works reasonable well", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat, cluster_idx_truth)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  dat_impute <- .impute_matrix(dat_NA, ncol(b_mat), test_prop)
  cluster_idx <- .spectral_clustering(dat_impute, ncol(b_mat))
  res <- .form_prediction_sbm(dat_impute, cluster_idx)
  
  b_mat2 <- matrix(0.7, 3, 3)
  b_mat2 <- b_mat2 - 0.5*diag(3)
  dat2 <- generate_sbm(b_mat2, cluster_idx)
  dat_NA2 <- .remove_entries(dat2, test_prop)
  dat_impute2 <- .impute_matrix(dat_NA2, ncol(b_mat2), test_prop)
  cluster_idx2 <- .spectral_clustering(dat_impute2, ncol(b_mat2))
  res2 <- .form_prediction_sbm(dat_impute2, cluster_idx2)
  
  expect_true(sum(abs(res - b_mat)) <= sum(abs(res2 - b_mat)))
})

################

## .sbm_projection is correct

test_that(".sbm_projection works", {
  set.seed(10)
  b_mat_truth <- 0.5*diag(3)
  b_mat_truth <- b_mat_truth + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  dat_impute <- .impute_matrix(dat_NA, ncol(b_mat_truth), test_prop)
  
  res <- .sbm_projection(dat_impute, ncol(b_mat_truth))
  
  expect_true(all(dim(res) == dim(dat)))
})

test_that(".sbm_projection works reasonably well", {
  set.seed(10)
  b_mat_truth <- 0.5*diag(3)
  b_mat_truth <- b_mat_truth + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  test_prop <- 0.1
  dat_NA <- .remove_entries(dat, test_prop)
  dat_impute <- .impute_matrix(dat_NA, ncol(b_mat_truth), test_prop)
  res <- .sbm_projection(dat_impute, ncol(b_mat_truth))
  
  b_mat_truth2 <- matrix(0.7, 3, 3)
  b_mat_truth2 <- b_mat_truth2 - 0.5*diag(3)
  dat2 <- generate_sbm(b_mat_truth2, cluster_idx)
  dat_NA2 <- .remove_entries(dat2, test_prop)
  dat_impute2 <- .impute_matrix(dat_NA2, ncol(b_mat_truth2), test_prop)
  res2 <- .sbm_projection(dat_impute2, ncol(b_mat_truth2))
  
  cluster_mat <- .convert_idx_to_mat(cluster_idx_truth)
  pop_mat <- cluster_mat %*% b_mat_truth %*% t(cluster_mat)
  
  expect_true(sum(abs(res - pop_mat)) <= sum(abs(res2 - pop_mat)))
})

######################

## network_cv_sbm_sample_split is correct

test_that("network_cv_sbm_sample_split works", {
  set.seed(10)
  b_mat_truth <- 0.5*diag(3)
  b_mat_truth <- b_mat_truth + 0.2
  cluster_idx_truth <- rep(1:3, each = 50)
  dat <- generate_sbm(b_mat_truth, cluster_idx_truth)
  
  res <- network_cv_sbm_sample_split(dat, k_vec = c(1:5), test_prop = 0.1)
  
  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(sort(names(res)) == sort(c("k", "err_vec", "err_mat"))))
  expect_true(length(res$err_vec) == 5)
  expect_true(ncol(res$err_mat) == length(res$err_vec))
})
