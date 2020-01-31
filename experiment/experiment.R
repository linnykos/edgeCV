rm(list = ls())
set.seed(20)
b_mat <- 0.5*diag(3)
b_mat <- b_mat + 0.2

b_tensor <- array(NA, c(10, 3, 3))
for(i in 1:dim(b_tensor)[1]){
  b_tensor[i,,] <- b_mat
}

cluster_idx <- rep(1:3, each = 20)

dat <- generate_tensor(b_tensor, cluster_idx, 1)
# res <- edge_cv_sbm_tensor(dat, 1:5, verbose = F)

k_vec = 1:5
nfolds = 5
tol = 1e-6
verbose = T

p <- dim(dat)[1]
n <- dim(dat)[2]

# assign each node triplet to a fold, stratified by first dimension
combn_mat <- utils::combn(n, 2)
idx_mat <- apply(combn_mat, 2, function(vec){
  c(.convert_pair_to_idx(vec, n), .convert_pair_to_idx(rev(vec), n))
})

fold_list <- lapply(1:p, function(j){
  fold_id <- rep(NA, ncol(idx_mat))
  for(i in 1:(nfolds-1)){
    fold_id[sample(which(is.na(fold_id)), round(ncol(idx_mat)/nfolds))] <- i
  }
  fold_id[which(is.na(fold_id))] <- nfolds
  
  fold_id
})


