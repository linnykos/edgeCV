rm(list=ls())
set.seed(10)
b_mat_truth <- 0.5*diag(3)
b_mat_truth <- b_mat_truth + 0.2
cluster_idx_truth <- rep(1:3, each = 50)
dat <- generate_sbm(b_mat_truth, cluster_idx_truth)

k_vec = c(1:5)
nfold = 5
tol = 1e-6
verbose = T

stopifnot(nrow(dat) == ncol(dat), sum(abs(dat - t(dat))) <= tol, 
          k_vec >= 0)

# assign each node pair to a fold
n <- nrow(dat)
combn_mat <- utils::combn(n, 2)
idx_mat <- apply(combn_mat, 2, function(vec){
  c(.convert_pair_to_idx(vec, n), .convert_pair_to_idx(rev(vec), n))
})
fold_id <- rep(NA, ncol(idx_mat))
for(i in 1:(nfold-1)){
  fold_id[sample(which(is.na(fold_id)), round(ncol(idx_mat)/nfold))] <- i
}
fold_id[which(is.na(fold_id))] <- nfold

err_mat_list <- lapply(1:nfold, function(fold){
  if(verbose) print(paste0("On fold ", fold))
  
  test_idx <- as.numeric(idx_mat[,which(fold_id == fold)])
  dat_NA <- dat
  dat_NA[test_idx] <- NA
  
  # generate the error matrix
  err_mat <- matrix(NA, nrow = sum(is.na(dat_NA)), ncol = length(k_vec))
  colnames(err_mat) <- k_vec
  
  # impute and compute errors
  for(i in k_vec){
    dat_impute <- .impute_matrix(dat_NA, k_vec[i], 1/nfold)
    dat_impute <- .sbm_projection(dat_impute, k_vec[i])
    err_mat[,i] <- (dat_impute[test_idx] - dat[test_idx])^2
  }
  
  err_mat
})

err_mat <- do.call(rbind, err_mat_list)

# compute the best model
err_vec <- colMeans(err_mat)
k <- k_vec[which.min(err_vec)]
