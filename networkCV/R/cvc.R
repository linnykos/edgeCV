cvc_sbm_sample_split <- function(err_mat, trials, alpha, verbose = T){
  stopifnot(length(colnames(err_mat)) == ncol(err_mat))
  
  n <- nrow(err_mat)
  k <- ncol(err_mat)
  err_mat2 <- .clean_err_mat(err_mat)
  
  mu_vec <- colMeans(err_mat2)
  sd_vec <- apply(err_mat2, 2, sd)
  test_vec <- .extract_max(sqrt(n)*mu_vec/sd_vec, k)
  
  boot_mat <- sapply(1:trials, function(b){
    if(verbose && b %% floor(trials/10) == 0) cat('*')
    tmp_mat <- .cvc_bootstrap_trial(err_mat2, mu_vec, sd_vec, k)
    .extract_max((1/sqrt(n)) * colSums(tmp_mat), k)
  })
  
  p_vec <- sapply(1:k, function(x){
    length(which(boot_mat[x,] <= test_vec[x]))/ncol(boot_mat)
  })
  
  model_selected <- colnames(err_mat)[which(p_vec <= alpha)]
  
  list(model_selected = model_selected, p_vec = p_vec)
}

###########

.clean_err_mat <- function(err_mat){
  k <- ncol(err_mat)
  err_mat2 <- do.call(cbind, lapply(1:k, function(x){
    sapply(c(1:k)[-x], function(y){
      err_mat[,x] - err_mat[,y]
    })
  }))
  
  colnames(err_mat2) <- unlist(lapply(1:k, function(x){sapply(c(1:k)[-x], function(y){paste0(x,"-",y)})}))
  err_mat2
}

.cvc_bootstrap_trial <- function(err_mat2, mu_vec, sd_vec, k){
  n <- nrow(err_mat2)
  g_vec <- stats::rnorm(n)
  ((diag(g_vec) %*% err_mat2) - (g_vec %*% t(mu_vec)))/(rep(1, n) %*% t(sd_vec))
}

.extract_max <- function(vec, k){
  stopifnot(length(vec) == k * (k-1))
  sapply(1:k, function(x){
    max(vec[((x-1)*(k-1)+1) : (x*(k-1))])
  })
}

