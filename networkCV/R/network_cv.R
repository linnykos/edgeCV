network_cv_sbm_sample_split <- function(dat, k_vec, test_prop, tol = 1e-6){
  stopifnot(nrow(dat) == ncol(dat), sum(abs(dat - t(dat))) <= tol, test_prop > 0, test_prop < 1,
            k_vec >= 0)

  # generate missing values
  dat_NA <- .remove_entries(dat, test_prop)
  test_idx <- which(is.na(dat_NA))
  
  # generate the error matrix
  err_mat <- matrix(NA, nrow = sum(is.na(dat_NA)), ncol = length(k_vec))
  colnames(err_mat) <- k_vec
  
  # impute and compute errors
  for(i in k_vec){
    dat_impute <- .impute_matrix(dat_NA, k_vec[i], test_prop)
    dat_impute <- .sbm_projection(dat_impute, k_vec[i])
    err_mat[,i] <- (dat_impute[test_idx] - dat[test_idx])^2
  }
  
  # compute the best model
  err_vec <- colMeans(err_mat)
  k <- k_vec[which.min(err_vec)]
  
  # output
  list(k = k, err_vec = err_vec, err_mat = err_mat)
}

###

.remove_entries <- function(dat, test_prop){
  n <- nrow(dat)
  dat_NA <- dat
  
  combn_mat <- utils::combn(n, 2)
  idx <- sample(1:ncol(combn_mat), round(test_prop * ncol(combn_mat)))
  
  for(i in 1:length(idx)){
    dat_NA[combn_mat[1,idx[i]], combn_mat[2,idx[i]]] <- NA
    dat_NA[combn_mat[2,idx[i]], combn_mat[1,idx[i]]] <- NA
  }
  
  dat_NA
}

.impute_matrix <- function(dat_NA, K, test_prop){
  test_idx <- which(is.na(dat_NA))
  dat_tmp <- dat_NA
  if(length(test_idx) > 0) dat_tmp[test_idx] <- 0
  eigen_res <- mgcv::slanczos(dat_tmp/(1-test_prop), K)
  
  if(K == 1){
    diag_mat <- matrix(eigen_res$values, 1, 1)
  } else {
    diag_mat <- diag(eigen_res$values)
  }
  
  eigen_res$vectors %*% diag_mat %*% t(eigen_res$vectors)
}

.spectral_clustering <- function(dat_impute, K){
  eigenvectors <- mgcv::slanczos(dat_impute, K)$vectors
  
  stats::kmeans(eigenvectors, centers=K, nstart=20)$cluster
}

.form_prediction_sbm <- function(dat_impute, cluster_idx){
  diag(dat_impute) <- 0
  K <- max(cluster_idx)
  
  b_mat <- matrix(NA, K, K)
  for(i in 1:K){
    for(j in 1:i){
      idx_i <- which(cluster_idx == i)
      
      if(i != j){
        idx_j <- which(cluster_idx == j)
        b_mat[i,j] <- mean(dat_impute[idx_i, idx_j])
        b_mat[j,i] <- b_mat[i,j]
        
      } else {
        b_mat[i,i] <- sum(dat_impute[idx_i, idx_i])/(2*choose(length(idx_i), 2))
      }
    }
  }
  
  b_mat
}

.sbm_projection <- function(dat_impute, K){
  cluster_idx <- .spectral_clustering(dat_impute, K)
  b_mat <- .form_prediction_sbm(dat_impute, cluster_idx)
  
  sapply(1:ncol(dat_impute), function(x){
    b_mat[cluster_idx[x], cluster_idx]
  })
}

