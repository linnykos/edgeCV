tensor_clustering <- function(dat, K, reps = 10, maxit = 100, verbose = T, tol = 1e-6){
  stopifnot(class(dat) == "array", length(dim(dat)) == 3,
            dim(dat)[2] == dim(dat)[3])
  
  p <- dim(dat)[1]
  n <- dim(dat)[2]
  
  for(i in 1:p){
    stopifnot(sum(abs(dat[i,,] - t(dat[i,,]))) <= 1e-6)
  }
  
  ss_best <- Inf
  clustering_best <- NULL
  b_tensor_best <- NULL
  
  dat_vec <- t(matrix(dat, p*n, n))
  
  for(rep in 1:reps){
    clustering <- kmeans(dat_vec, K)$cluster 
    iter <- 0
    ss_prev <- Inf
    
    while (iter < maxit){
      if (verbose) cat('step = ', iter, '\n' )
      iter <- iter + 1
      
      # deal with cluster do not have any members
      counts <- table(clustering)
      if(any(counts == 0)) clustering <- .resolve_empty_cluster(clustering)
      
      b_tensor <- .form_prediction_tensor(dat, clustering)
      dist_mat <- .compute_distance_tensors(dat, b_tensor)
      ss_current <- .compute_ss(dist_mat)
      
      if (ss_current - ss_prev > 1e-6) break()
      
      clustering_prev <- clustering
      ss_prev <- ss_current
      b_tensor_prev <- b_tensor
      
      clustering <- .update_cluster_tensors(clustering, dist_mat)
    }
    
    if (ss_current < ss_best){
      ss_best <- ss_current
      clustering_best <- clustering
      b_tensor_best <- b_tensor
    }
  }
  
  list(clustering = clustering_best, b_tensor = b_tensor_best)
}

.resolve_empty_clustering <- function(vec){
  empties <- which(counts == 0)
  for (i in empties){
    from <- sample(which(counts > 1), 1,0)
    lonely <- sample(which(idx == from), 1, 0)
    idx[lonely] <- i
    counts[i] <- 1
    counts[from] <- counts[from]-1; 
  }   
}


.form_prediction_tensor <-function(A, idx, k){
  p <- dim(A)[1]
  m.1 <- dim(A)[2]
  idx.1 <- idx[1:m.1]
  center <- array(0, dim=c(p,k,k))
  for (i in 1:k){
    for (j in i:k){
      if ((sum(idx.1 == i) ==1) & (sum(idx==j)==1) ) {
        center[,i,j] <-0.5
        center[,j,i] <- center[,i,j]
      } else{
        center[,i,j] <- apply(A[,idx.1==i,idx==j],1,mean)
        center[,j,i] <- center[,i,j]
      }
    }
  }
  return(center)
}

.compute_distance_tensors <- function(dat, b_tensor){
  # A: p*m*m
  # center: p*m*k
  k <- dim(center)[3]
  m <- dim(A)[3]
  D <- matrix(0, m, k)
  for (j in 1:m){
    for (i in 1:k){
      D[j,i] <- norm( (A[,,j] - center[,,i]), 'F')^2
    }
  }
  return(D)
}

.compute_ss <- function(){
  totsumD <- sum(D[(idx-1)*m+(1:m)])
}


.update_cluster_tensors <- function(clustering, dist_mat){
  nidx <- apply(D, 1, which.min)
  moved <- which(nidx != previdx)
  # resolve tie in favor of not move
  if (length(moved)>0) {
    moved <- moved[D[(previdx[moved]-1)*m+moved] > min(D[moved,])] 
  } 
  if (length(moved) ==0){
    break
  }
  idx[moved] = nidx[moved]
}
