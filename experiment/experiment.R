rm(list=ls())
library(simulation)
library(networkCV)

set.seed(10)
trials <- 50
paramMat <- as.matrix(expand.grid(c(600), c(seq(0, 0.6, length.out = 13)),
                                  c(3:5),
                                  6 , 200, 5))
colnames(paramMat) <- c("n", "rho", "K", "num_model", "trials", "nfold")

#############

rule <- function(vec){
  b_mat <- matrix(0.2, vec["K"], vec["K"]) + 0.6*diag(vec["K"])
  n <- vec["n"]
  n_each <- round(n/vec["K"])
  cluster_idx <- rep(1:(vec["K"]-1), each = n_each)
  cluster_idx <- c(cluster_idx, rep(vec["K"], n - length(cluster_idx)))
  if(vec["rho"] >= 0){
    rho <- 1/(n^vec["rho"])
  } else {
    rho <- log(n)/n
  }
  dat <- networkCV::generate_sbm(b_mat, cluster_idx, rho)
  
  dat
}

set.seed(1)
vec = paramMat[1,]
dat = rule(vec)
k_vec = c(1:vec["num_model"])
nfolds = vec["nfold"]
verbose = F
tol = 1e-6

stopifnot(nrow(dat) == ncol(dat), sum(abs(dat - t(dat))) <= tol, 
          k_vec >= 0)

# assign each node pair to a fold
n <- nrow(dat)
combn_mat <- utils::combn(n, 2)
idx_mat <- apply(combn_mat, 2, function(vec){
  c(.convert_pair_to_idx(vec, n), .convert_pair_to_idx(rev(vec), n))
})
fold_id <- rep(NA, ncol(idx_mat))
for(i in 1:(nfolds-1)){
  fold_id[sample(which(is.na(fold_id)), round(ncol(idx_mat)/nfolds))] <- i
}
fold_id[which(is.na(fold_id))] <- nfolds

err_mat_list <- lapply(1:nfolds, function(fold){
  if(verbose) print(paste0("On fold ", fold))
  
  test_idx <- as.numeric(idx_mat[,which(fold_id == fold)])
  dat_NA <- dat
  dat_NA[test_idx] <- NA
  
  # generate the error matrix
  err_mat <- matrix(NA, nrow = sum(is.na(dat_NA)), ncol = length(k_vec))
  colnames(err_mat) <- k_vec
  
  # impute and compute errors
  for(i in k_vec){
    dat_impute <- .impute_matrix(dat_NA, k_vec[i], 1/nfolds)
    dat_impute <- .sbm_projection(dat_impute, k_vec[i], dat_org = NA)
    err_mat[,i] <- (dat_impute[test_idx] - dat[test_idx])^2
  }
  
  err_mat
})

err_mat <- do.call(rbind, err_mat_list)

# compute the best model
err_vec <- colMeans(err_mat)