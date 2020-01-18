context("Test generate_sbm")

## .convert_idx_to_mat is correct

test_that(".convert_idx_to_mat works", {
  cluster_idx <- rep(1:5, each = 5)
  res <- .convert_idx_to_mat(cluster_idx)
  
  expect_true(all(dim(res) == c(25,5)))
})

#####################

## generate_sbm is correct

test_that("generate_sbm works", {
  b_mat <- diag(3)
  cluster_idx <- rep(1:3, each = 5)
  res <- generate_sbm(b_mat, cluster_idx)
  
  expect_true(all(dim(res) == 15))
  expect_true(sum(abs(res - t(res))) <= 1e-6)
})

test_that("generate_sbm works for a random matrix", {
  trials <- 100
  
  b_mat1 <- 0.5*diag(3)
  b_mat1 <- b_mat1 + 0.2
  cluster_idx1 <- rep(1:3, each = 5)
  res1_list <- lapply(1:trials, function(x){generate_sbm(b_mat1, cluster_idx1)})
  res1 <- Reduce("+", res1_list) / length(res1_list)
  
  b_mat2 <- matrix(0.7, 3, 3)
  b_mat2 <- b_mat2 - 0.5*diag(3)
  cluster_idx2 <- rep(1:3, each = 5)
  res2_list <- lapply(1:trials, function(x){generate_sbm(b_mat2, cluster_idx2)})
  res2 <- Reduce("+", res2_list) / length(res2_list)
  
  cluster_mat1 <- .convert_idx_to_mat(cluster_idx1)
  pop_mat1 <- cluster_mat1 %*% b_mat1 %*% t(cluster_mat1)
  
  expect_true(sum(abs(pop_mat1 - res1)) <= sum(abs(pop_mat1 -res2)))
})
