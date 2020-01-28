tensor_clustering <- function(dat, K, reps = 10, maxit = 100, verbose = F, tol = 1e-6){
  stopifnot(class(dat) == "array", length(dim(dat)) == 3,
            dim(dat)[2] == dim(dat)[3])
  
  p <- dim(dat)[1]
  n <- dim(dat)[2]
  
  for(i in 1:p){
    stopifnot(sum(abs(dat[i,,] - t(dat[i,,]))) <= 1e-6)
  }
  
  if(K == 1){
    clustering_best <- rep(1, n)
    b_tensor_best <- array(mean(dat), c(1,1,1))
    
    return(list(clustering = clustering_best, b_tensor = b_tensor_best))
  }
  
  ss_best <- Inf
  clustering_best <- NULL
  b_tensor_best <- NULL
  
  dat_vec <- t(matrix(dat, p*n, n))
  
  for(rep in 1:reps){
    clustering <- stats::kmeans(dat_vec, K)$cluster 
    iter <- 0
    ss_prev <- Inf
    
    while(iter < maxit){
      if(verbose) cat('step = ', iter, '\n' )
      iter <- iter + 1
      
      # deal with cluster do not have any members
      clustering <- .resolve_empty_cluster(clustering, K)
      
      b_tensor <- .form_prediction_tensor(dat, clustering)
      dist_mat <- .compute_distance_tensors(dat, b_tensor, clustering)
      ss_current <- .compute_ss(clustering, dist_mat)
      
      if(ss_current - ss_prev > 1e-6) break()
      
      clustering_prev <- clustering
      ss_prev <- ss_current
      b_tensor_prev <- b_tensor
      
      clustering <- .update_cluster_tensors(dist_mat)
      if(all(clustering == clustering_prev)) break()
    }
    
    if(ss_current < ss_best){
      ss_best <- ss_current
      clustering_best <- clustering
      b_tensor_best <- b_tensor
    }
  }
  
  list(clustering = clustering_best, b_tensor = b_tensor_best)
}

.resolve_empty_cluster <- function(clustering, K){
  counts <- table(factor(clustering, levels = 1:K))
  tmp <- clustering
  
  if(any(counts == 0)) {
    empty_idx <- which(counts == 0)
    
    if(length(empty_idx) > 0){
      
      for (i in empty_idx){
        # move idx into the empty cluster
        from <- sample(which(counts > 1), 1)
        idx <- sample(which(clustering == from), 1)
        tmp[idx] <- i
        counts[i] <- 1
        counts[from] <- counts[from]-1
      }
    }
  }
    
  tmp
}


.form_prediction_tensor <-function(dat, clustering){
  K <- max(clustering)
  
  p <- dim(dat)[1]
  n <- dim(dat)[2]
  
  b_tensor <- array(NA, dim=c(p,K,K))
  
  for (i in 1:K){
    for (j in i:K){
      # only one node in both clusters
      if ((sum(clustering == i) == 1) & (sum(clustering == j) == 1)) {
        b_tensor[,i,j] <- 0.5
        b_tensor[,j,i] <- b_tensor[,i,j]
        
      } else{
        b_tensor[,i,j] <- apply(dat[,clustering==i, clustering==j], 1, mean)
        b_tensor[,j,i] <- b_tensor[,i,j]
      }
    }
  }

  b_tensor
}

.l2norm <- function(x){sqrt(sum(x^2))}

.compute_distance_tensors <- function(dat, b_tensor, clustering){
  stopifnot(length(clustering) == dim(dat)[2], all(clustering <= dim(b_tensor)[2]))
  # dat: p*n*n
  # b_tensor: p*k*k
  
  K <- dim(b_tensor)[3]
  n <- dim(dat)[3]
  tmp <- b_tensor[,clustering,]
  dist_mat <- matrix(0, n, K)
  
  for (j in 1:n){
    for (i in 1:K){
      dist_mat[j,i] <- .l2norm(dat[,,j] - tmp[,,i])
    }
  }
  
  dist_mat
}

# uses fancy indexing to grab the value of dist_mat corresponding to each point
.compute_ss <- function(clustering, dist_mat){
  n <- nrow(dist_mat)
  
  sum(dist_mat[(clustering-1)*n+(1:n)])
}

.update_cluster_tensors <- function(dist_mat){
  apply(dist_mat, 1, which.min)
}
