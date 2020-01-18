generate_sbm <- function(b_mat, cluster_idx){
  stopifnot(max(cluster_idx) == nrow(b_mat), sum(abs(b_mat - t(b_mat))) <= 1e-6)
  
  n <- length(cluster_idx)
  cluster_mat <- .convert_idx_to_mat(cluster_idx)
  pop_mat <- cluster_mat %*% b_mat %*% t(cluster_mat)
  
  adj_mat <- matrix(0, n, n)
  lower_tri_idx <- which(lower.tri(adj_mat, diag = F))
  adj_mat[lower_tri_idx] <- stats::rbinom(length(lower_tri_idx), size = 1, prob = pop_mat[lower_tri_idx])
  adj_mat <- t(adj_mat) + adj_mat
  
  adj_mat
}

.convert_idx_to_mat <- function(cluster_idx){
  k <- max(cluster_idx)
  n <- length(cluster_idx)
  
  sapply(1:k, function(i){
    tmp <- rep(0, n)
    tmp[cluster_idx == i] <- 1
    tmp
  })
}