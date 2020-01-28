rm(list=ls())
set.seed(20)
b_mat <- 0.5*diag(3)
b_mat <- b_mat + 0.2

b_tensor <- array(NA, c(10, 3, 3))
for(i in 1:dim(b_tensor)[1]){
  b_tensor[i,,] <- b_mat
}

cluster_idx <- rep(1:3, each = 5)

dat <- generate_tensor(b_tensor, cluster_idx, 1)
k_vec <- 1:5
trials = 5
test_prop = 0.1
tol = 1e-6
verbose = T

p <- dim(dat)[1]
n <- dim(dat)[2]

test_idx <- .generate_tensor_indices(p, n, round(test_prop*p*n*(n-1)/2))

dat_NA <- dat
dat_NA[test_idx] <- NA

# generate the error matrix
err_mat <- matrix(NA, nrow = sum(is.na(dat_NA)), ncol = length(k_vec))
colnames(err_mat) <- k_vec

# impute and compute errors
for(i in k_vec){
  dat_impute <- .impute_tensor(dat_NA, k_vec[i], test_prop)
  err_mat[,i] <- (dat_impute[test_idx] - dat[test_idx])^2
}

