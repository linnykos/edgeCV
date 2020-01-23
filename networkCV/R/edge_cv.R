edge_cv_sbm_sample_split <- function(dat, k_vec, test_prop, tol = 1e-6){
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

edge_cv_sbm <- function(dat, k_vec, nfold = 5, tol = 1e-6, verbose = T){
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
  
  # output
  list(k = k, err_vec = err_vec, err_mat_list = err_mat_list)
}

###

.convert_pair_to_idx <- function(vec, n){
  (vec[1]-1)*n+vec[2]
}

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

