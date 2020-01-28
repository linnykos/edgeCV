context("Test tensor clustering")

## .resolve_empty_cluster is correct

test_that(".resolve_empty_cluster works", {
  set.seed(10)
  clustering <- sample(1:5, 100, replace = T)
  res <- .resolve_empty_cluster(clustering, 5)
  
  expect_true(all(clustering == res))
})

test_that(".resolve_empty_cluster actually does something", {
  set.seed(10)
  clustering <- sample(1:5, 100, replace = T)
  clustering[clustering == 2] <- 1
  res <- .resolve_empty_cluster(clustering, 5)
  
  expect_true(sum(clustering != res) == 1)
  expect_true(length(which(res == 2)) == 1)
})

test_that(".resolve_empty_cluster works when the missing cluster is the last one", {
  set.seed(10)
  clustering <- sample(1:4, 100, replace = T)
  res <- .resolve_empty_cluster(clustering, 5)
  
  expect_true(sum(clustering != res) == 1)
  expect_true(length(which(res == 5)) == 1)
})

##########################

## .form_prediction_tensor is correct

test_that(".form_prediction_tensor works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  clustering <- rep(1:3, each = 5)
  
  dat <- generate_tensor(b_tensor, clustering, 1)
  
  res <- .form_prediction_tensor(dat, clustering)
  
  expect_true(all(dim(res) == dim(b_tensor)))
})

###################

## .compute_distance_tensors is correct

test_that(".compute_distance_tensors works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  clustering <- rep(1:3, each = 5)
  
  dat <- generate_tensor(b_tensor, clustering, 1)
  
  b_est <- .form_prediction_tensor(dat, clustering)
  res <- .compute_distance_tensors(dat, b_est, clustering)
  
  expect_true(all(dim(res) == c(dim(dat)[2], dim(b_tensor)[2])))
})

##################

## .compute_ss is correct

test_that(".compute_ss works", {
  set.seed(10)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  clustering <- rep(1:3, each = 5)
  
  dat <- generate_tensor(b_tensor, clustering, 1)
  
  b_est <- .form_prediction_tensor(dat, clustering)
  dist_mat <- .compute_distance_tensors(dat, b_est, clustering)
  
  res <- .compute_ss(clustering, dist_mat)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".compute_ss gives the correct calculation", {
  set.seed(20)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  clustering <- rep(1:3, each = 5)
  
  dat <- generate_tensor(b_tensor, clustering, 1)
  
  b_est <- .form_prediction_tensor(dat, clustering)
  dist_mat <- .compute_distance_tensors(dat, b_est, clustering)
  
  res <- .compute_ss(clustering, dist_mat)
  
  res2 <- 0
  for(i in 1:nrow(dist_mat)){
    res2 <- res2 + dist_mat[i,clustering[i]]
  }
  
  expect_true(abs(res - res2) <= 1e-6)
})

#####################

## tensor_clustering is correct

test_that("tensor_clustering works", {
  set.seed(20)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  clustering <- rep(1:3, each = 100)
  
  dat <- generate_tensor(b_tensor, clustering, 1)
  
  res <- tensor_clustering(dat, 3, verbose = F)
  
  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("clustering", "b_tensor")))
})

test_that("tensor_clustering works for K=1", {
  set.seed(20)
  b_mat <- 0.5*diag(3)
  b_mat <- b_mat + 0.2
  
  b_tensor <- array(NA, c(10, 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  clustering <- rep(1:3, each = 100)
  
  dat <- generate_tensor(b_tensor, clustering, 1)
  
  res <- tensor_clustering(dat, 1, verbose = F)
  
  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("clustering", "b_tensor")))
})
