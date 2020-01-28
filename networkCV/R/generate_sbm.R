#' Generate SBM
#'
#' @param b_mat connectivity matrix (symmetric), all values between \code{0} and \code{1} 
#' @param cluster_idx vector of length equal to number of nodes, where each node is 
#' assigned to a positive integer from \code{1} to \code{max(cluster_idx)}
#' @param rho sparsity index, between \code{0} and \code{1} 
#'
#' @return symmetric \code{0}-\code{1} matrix
#' @export
generate_sbm <- function(b_mat, cluster_idx, rho = 1){
  stopifnot(max(cluster_idx) == nrow(b_mat), sum(abs(b_mat - t(b_mat))) <= 1e-6,
            rho >= 0, rho <= 1)
  
  n <- length(cluster_idx)
  cluster_mat <- .convert_idx_to_mat(cluster_idx)
  pop_mat <- rho*(cluster_mat %*% b_mat %*% t(cluster_mat))
  
  adj_mat <- matrix(0, n, n)
  lower_tri_idx <- which(lower.tri(adj_mat, diag = F))
  adj_mat[lower_tri_idx] <- stats::rbinom(length(lower_tri_idx), size = 1, prob = pop_mat[lower_tri_idx])
  adj_mat <- t(adj_mat) + adj_mat
  
  adj_mat
}

#' Generate tensor SBM
#'
#' @param b_tensor connectivity tensor as an \code{array}, organized so the dimension is \code{p} by \code{n} by \code{n},
#' representing \code{p} symmetric matrices all values between \code{0} and \code{1} 
#' @param cluster_idx vector of length equal to the number of nodes \code{n}, where each node is 
#' assigned to a positive integer from \code{1} to \code{max(cluster_idx)}
#' @param rho  sparsity index, between \code{0} and \code{1} 
#'
#' @return a \code{array}
#' @export
generate_tensor <- function(b_tensor, cluster_idx, rho) {
  stopifnot(class(b_tensor) == "array", length(dim(b_tensor)) == 3,
            dim(b_tensor)[2] == dim(b_tensor)[3])
  
  k <- dim(b_tensor)[2]
  n <- length(cluster_idx)
  p <- dim(b_tensor)[1]
  
  adj_array <- array(0, dim = c(p, n, n))
  
  for (i in 1:p){ 
    adj_array[i,,] <- generate_sbm(b_tensor[i,,], cluster_idx, rho = rho)
  }
  
  
  adj_array
}

##################3

.convert_idx_to_mat <- function(cluster_idx){
  k <- max(cluster_idx)
  n <- length(cluster_idx)
  
  sapply(1:k, function(i){
    tmp <- rep(0, n)
    tmp[cluster_idx == i] <- 1
    tmp
  })
}