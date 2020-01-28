#' Edge CV for tensors (cross validation)
#'
#' @param dat \code{p} by \code{n} by \code{n} symmetric adjacency matrix
#' @param k_vec vector of positive integers, one per model, of how many blocks in said model
#' @param trials positive integer
#' @param test_prop number between 0 and 1
#' @param tol small positive number
#' @param verbose boolean
#'
#' @return list
#' @export
edge_cv_sbm_tensor <- function(dat, k_vec, trials = 5, test_prop = 0.1, tol = 1e-6, verbose = T){
  p <- dim(dat)[1]
  n <- dim(dat)[2]
  
  err_mat_list <- lapply(1:trials, function(trial){
    if(verbose) print(paste0("On trial ", trial))
    
    test_idx <- .generate_tensor_indices(p, n, round(test_prop*p*n*(n+1)/2))
    
    dat_NA <- dat
    dat_NA[test_idx] <- NA
    
    # generate the error matrix
    err_mat <- matrix(NA, nrow = length(test_idx), ncol = length(k_vec))
    colnames(err_mat) <- k_vec
    
    # impute and compute errors
    for(i in k_vec){
      dat_impute <- .impute_tensor(dat_NA, k_vec[i], test_prop)
      err_mat[,i] <- (dat_impute[test_idx] - dat[test_idx])^2
    }
    
    err_mat
  })
  
  err_mat <- do.call(rbind, err_mat_list)
  
  # compute the best model
  err_vec <- colMeans(err_mat)
  k <- k_vec[which.min(err_vec)]
  
  # output
  list(k = k, err_vec = err_vec, err_mat_list = err_mat_list)
}

#############

.generate_tensor_indices <- function(p, n, num){
  combn_mat <- utils::combn(n, 2)
  
  tmp_mat <- matrix(NA, ncol = 3, nrow = num)
  for(i in 1:num){
    tmp_mat[i,1] <- sample(1:p, 1)
    tmp_mat[i,2:3] <- sample(combn_mat[,sample(1:ncol(combn_mat),1)])
  }
  
  vec1 <- apply(tmp_mat, 1, .convert_triplet_to_idx, p = p, n = n)
  vec2 <- apply(tmp_mat[,c(1,3,2)], 1, .convert_triplet_to_idx, p = p, n = n)
  
  c(vec1, vec2)
}

.convert_triplet_to_idx <- function(vec, p, n){
  (vec[3]-1)*(p*n) + (vec[2]-1)*p + vec[1]
}

.impute_tensor <- function(dat_NA, K, test_prop){
  test_idx <- which(is.na(dat_NA))
  dat_tmp <- dat_NA
  if(length(test_idx) > 0) dat_tmp[test_idx] <- 0
  res <- tensor_clustering(dat_tmp/(1-test_prop), K)
  
  b_mat <- apply(res$b_tensor, c(2,3), mean)
  cluster_mat <- .convert_idx_to_mat(res$clustering)
  
  dat_impute <- array(NA, dim(dat_NA))
  for(i in 1:dim(dat_NA)[1]){
    dat_impute[i,,] <- cluster_mat %*% b_mat %*% t(cluster_mat)
  }
  
  dat_impute
}
