NCV <- function(A, K.vec = 1:6, n.fold = 5, n.dim = NULL) {
  # type = type[1]
  # method = method[1]
  n <- nrow(A)
  sample.ind <- sample.int(n)
  fold.size <- floor(n / n.fold)
  fold.ind <- list()
  for (i in 1:(n.fold - 1)) {
    fold.ind[[i]] <- sort(sample.ind[((i-1)*fold.size+1):(i*fold.size)])
  }
  fold.ind[[n.fold]] <- sort(sample.ind[((n.fold - 1) * fold.size + 1) : n])
  
  #  n.err <- sum(sapply(fold.ind,function(x){length(x)*(length(x)-1)/2}))
  n.K <- length(K.vec)
  err.mat <- matrix(0, nrow = 0, ncol = n.K)
  
  for (i in 1:n.fold) {
    tmp.ind <- fold.ind[[i]]
    A.1 <- A[-tmp.ind, -tmp.ind]
    A.12 <- A[-tmp.ind, tmp.ind]
    A.2 <- A[tmp.ind, tmp.ind]
    err.mat.tmp <- lapply(1:n.K, function(j){Validate.CVC.Rec(A.1, K.vec[j], 
                                                              A.12, A.2, n.dim = ifelse(is.null(n.dim), K.vec[j], n.dim))})
    err.mat.tmp <- do.call(cbind, err.mat.tmp)
    err.mat <- rbind(err.mat, err.mat.tmp)
  }
  
  # compute the best model
  err.mat <- err.mat^2
  err_vec <- colMeans(err.mat^2)
  k <- K.vec[which.min(err_vec)]
  
  # output
  list(k = k, err_vec = err_vec, err_mat = err.mat)
}

#########################


Validate.CVC.Rec <- function(A.1, K, A.12, A.2, n.dim = K){
  n.1 <- nrow(A.1)
  n.2 <- nrow(A.2)
  n <- n.1 + n.2
  A.rec <- cbind(A.1, A.12)
  
  if (K == 1) {
    #p <- (sum(A.1) + 2) / (n.1^2 - n.1 + 2)
    p <- (sum(A.1) / 2 + sum(A.12) + 1) / ((n.1^2 - n.1) / 2 + n.1 * n.2 + 1)
    pred.err.mat <- A.2 - p
    err <- as.vector(pred.err.mat[upper.tri(pred.err.mat)])
  } else {
    clusters <- SpecClust.Rec(A = t(A.rec), K = K, n.dim = n.dim)
    cluster.2 <- clusters[(n.1 + 1):n]  
    # if (length(unique(cluster.2)) < K) {
    #   err = Inf
    #   return(err)
    # } else {
    
    Theta.2 <- ClustVec2Mat(cluster.2, K)   
    B.rec <- EmpB.Rec(A.1, A.12, clusters)      
    P.2 <- Theta.2 %*% B.rec %*% t(Theta.2)
    pred.err.mat <- A.2 - P.2
    err <- as.vector(pred.err.mat[upper.tri(pred.err.mat)])
    # return(err)
    #}
  }
  return(err)
}


SpecClust.Rec <- function(A, K, sphere = F, n.dim = K) {
  # spectral clustering on rectangular
  n <- nrow(A)
  if (K == 1) {
    return(rep(1, n))
  }
  U <- svd(A, nv = 0, nu = n.dim)$u
  #U <- slanczos(A%*%t(A), n.dim)$vectors
  #U <- propack.svd(A, n.dim)$u
  #U <- irlba(A, nu = n.dim, nv = 0)$u
  if (sphere) {
    for (i in 1:n) {
      U[i,] = U[i,] / sqrt(sum(U[i,]^2))
    }
  }
  return(kmeans(U, K, nstart = 20)$cluster)
  #  return(clusters(kcca(U, K)))
}


ClustVec2Mat <- function(clust.vec, K) {
  # convert a membership vector to a membership matrix, the number of different values in
  # clust.vec must be at most K
  clust.mat <- matrix(0, ncol = K, nrow = length(clust.vec))
  for (i in 1:K) {
    clust.mat[clust.vec==i,i] <- 1
  }
  return(clust.mat)
}

EmpB.Rec <- function(A.1, A.12, clusters) {
  # empirical connectivity matrix
  n.1 <- nrow(A.12)
  n.2 <- ncol(A.12)
  n <- n.1 + n.2
  values <- sort(unique(clusters))
  K <- length(values)
  B <- matrix(0, nrow = K, ncol = K)
  clust.ind.1 <- list()
  clust.size.1 <- rep(0, K)
  clust.ind.2 <- list()
  clust.size.2 <- rep(0, K)
  for (i in 1:K) {
    clust.ind.1[[i]] <- which(clusters[1:n.1] == i)
    clust.ind.2[[i]] <- which(clusters[(n.1 + 1):n] == i)
    clust.size.1[i] <- length(clust.ind.1[[i]])
    clust.size.2[i] <- length(clust.ind.2[[i]])
  }
  for (i in 1:K) {
    e.tmp <- sum(A.1[clust.ind.1[[i]], clust.ind.1[[i]]]) / 2 +
      sum(A.12[clust.ind.1[[i]], clust.ind.2[[i]]])
    n.tmp <- (clust.size.1[i]^2 - clust.size.1[i]) / 2 +
      clust.size.1[i] * clust.size.2[i]
    B[i,i] <- (e.tmp + 1) / (n.tmp + 1)
    if (i < K) {
      for (j in (i+1):K) {
        e.tmp <- sum(A.1[clust.ind.1[[i]], clust.ind.1[[j]]]) +
          sum(A.12[clust.ind.1[[i]], clust.ind.2[[j]]])
        n.tmp <- clust.size.1[i] * (clust.size.1[j] + clust.size.2[j])
        B[i,j] <- (e.tmp + 1) / (n.tmp + 1)
        B[j,i] <- B[i,j]
      }
    }
  }
  return(B)
}


