context("Test edge CV for tensors")

## .generate_tensor_indices is correct

test_that(".generate_tensor_indices works", {
  set.seed(10)
  res <- .generate_tensor_indices(10, 5, 50)
  
  expect_true(all(res <= 5*10*10))
})

test_that(".generate_tensor_indices keeps valid indices", {
  set.seed(20)
  test_prop <- 0.1
  p <- 10
  n <- 15
  res <- .generate_tensor_indices(p, n, round(test_prop*p*n*(n-1)/2))
  
  expect_true(all(res <= 10*15*15))
})

#################

## .convert_triplet_to_idx is correct

test_that(".convert_triplet_to_idx works", {
  res <- .convert_triplet_to_idx(c(2,4,3), 5, 4)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".convert_triplet_to_idx is correct", {
  mat <- array(c(1:c(5*7*7)), c(5,7,7))
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    vec <- c(sample(1:5, 1), sample(1:7, 2, replace = T))
    res <- .convert_triplet_to_idx(vec, 5, 7)
    
    res == mat[vec[1],vec[2],vec[3]]
  })
  
  expect_true(all(bool_vec))
})

###################

## .impute_tensor is correct

test_that(".impute_tensor works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  cluster_idx <- rep(1:3, each = 5)
  
  dat <- generate_tensor(b_tensor, cluster_idx, 1)
  p <- dim(dat)[1]; n <- dim(dat)[2]
  test_prop <- 0.1
  
  test_idx <- .generate_tensor_indices(p, n, round(test_prop*p*n*(n+1)/2))
  
  dat_NA <- dat
  dat_NA[test_idx] <- NA
  
  res <- .impute_tensor(dat_NA, 3, test_prop)
  
  expect_true(all(dim(res) == dim(dat_NA)))
  expect_true(sum(is.na(res)) == 0)
})

test_that(".impute_tensor works for K=1", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  cluster_idx <- rep(1:3, each = 5)
  
  dat <- generate_tensor(b_tensor, cluster_idx, 1)
  p <- dim(dat)[1]; n <- dim(dat)[2]
  test_prop <- 0.1
  
  test_idx <- .generate_tensor_indices(p, n, round(test_prop*p*n*(n+1)/2))
  
  dat_NA <- dat
  dat_NA[test_idx] <- NA
  
  res <- .impute_tensor(dat_NA, 1, test_prop)
  
  expect_true(all(dim(res) == dim(dat_NA)))
  expect_true(sum(is.na(res)) == 0)
})

######################

## edge_cv_sbm_tensor is correct

test_that("edge_cv_sbm_tensor works", {
  set.seed(20)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  cluster_idx <- rep(1:3, each = 50)
  
  dat <- generate_tensor(b_tensor, cluster_idx, 1)
  
  res <- edge_cv_sbm_tensor(dat, 1:5, verbose = F)
  
  expect_true(is.list(res))
  expect_true(length(res) == 3)
})