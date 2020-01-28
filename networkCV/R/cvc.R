#' CVC for SBM (with cross validation)
#'
#' @param err_mat_list error matrix list outputted by \code{edge_cv_sbm$err_mat_list}
#' @param trials number of trials
#' @param alpha significant level
#' @param verbose boolean
#' @param ncores integer or \code{NA}
#'
#' @return list
#' @export
cvc_sbm <- function(err_mat_list, trials, alpha, verbose = T, ncores = NA){
  if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)
  
  for(i in 1:length(err_mat_list)){
    stopifnot(length(colnames(err_mat_list[[i]])) == ncol(err_mat_list[[i]]))
  }
  
  n_vec <- sapply(err_mat_list, function(err_mat){nrow(err_mat)})
  k_vec <- sapply(err_mat_list, function(err_mat){ncol(err_mat)})
  stopifnot(length(unique(k_vec)) == 1)
  k <- unique(k_vec)
  n <- sum(sapply(err_mat_list, nrow))
  
  err_mat_list2 <- lapply(err_mat_list, function(err_mat){
    .clean_err_mat(err_mat)
  })
  
  mu_vec_list <- lapply(err_mat_list2, function(err_mat2){colMeans(err_mat2)})
  
  err_mat_list2_recentered <- .recenter_list(err_mat_list2, mu_vec_list)
  
  mu_vec <- colMeans(do.call(rbind, mu_vec_list))
  sd_vec <- apply(do.call(rbind, err_mat_list2_recentered), 2, stats::sd)
  
  # compute test statistics
  test_vec <- .extract_max(sqrt(n)*mu_vec/sd_vec, k)
  
  # compute null distribution
  func2 <- function(b){
    if(verbose && b %% floor(trials/10) == 0) cat('*')
    
    set.seed(sum(sd_vec)*100+10*b)
    tmp_mat <- .cvc_bootstrap_trial_cv(err_mat_list2, mu_vec_list, sd_vec)
    tmp <- .extract_max((1/sqrt(n)) * colSums(tmp_mat), k)
    
    if(verbose) print(paste0("Trial ", b, ": Values ", paste0(round(tmp,2), collapse = ", ")))
    tmp
  }
  
  if(!is.na(ncores)){
    b <- 0 #debugging purposes
    boot_mat <- foreach::"%dopar%"(foreach::foreach(b = 1:trials), func2(b))

    if(verbose) print("Done bootstrapping")
    boot_mat <- do.call(cbind, boot_mat)
  } else {
    boot_mat <- sapply(1:trials, func2)
  }
  
  # SUSPECTED BUG HERE
  if(verbose) print("Computing p-values")
  # compute p-value
  p_vec <- sapply(1:k, function(x){
    length(which(boot_mat[x,] <= test_vec[x]))/ncol(boot_mat)
  })
  
  model_selected <- colnames(err_mat_list[[1]])[which(p_vec <= alpha)]
  
  list(model_selected = model_selected, p_vec = p_vec)
}

#' CVC for SBM (with sample splitting)
#'
#' @param err_mat error matrix outputted by \code{edge_cv_sbm_sample_split$err_mat}
#' @param trials number of trials
#' @param alpha significant level
#' @param verbose boolean
#' @param ncores integer or \code{NA}
#'
#' @return list 
#' @export
cvc_sbm_sample_split <- function(err_mat, trials, alpha, verbose = T, ncores = NA){
  if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)
  
  stopifnot(length(colnames(err_mat)) == ncol(err_mat))
  
  n <- nrow(err_mat)
  k <- ncol(err_mat)
  err_mat2 <- .clean_err_mat(err_mat)
  
  mu_vec <- colMeans(err_mat2)
  sd_vec <- apply(err_mat2, 2, stats::sd)
  test_vec <- .extract_max(sqrt(n)*mu_vec/sd_vec, k)
  
  func <- function(b){
    if(verbose && b %% floor(trials/10) == 0) cat('*')
    
    set.seed(sum(sd_vec)*100+10*b)
    tmp_mat <- .cvc_bootstrap_trial(err_mat2, mu_vec, sd_vec)
    tmp <- .extract_max((1/sqrt(n)) * colSums(tmp_mat), k)
    
    if(verbose) print(paste0("Trial ", b, ": Values ", paste0(round(tmp,2), collapse = ", ")))
    tmp
  }
  
  if(!is.na(ncores)){
    b <- 0 #debugging purposes
    boot_mat <- foreach::"%dopar%"(foreach::foreach(b = 1:trials), func(b))
    boot_mat <- do.call(cbind, boot_mat)
  } else {
    boot_mat <- sapply(1:trials, func)
  }
  
  
  p_vec <- sapply(1:k, function(x){
    length(which(boot_mat[x,] <= test_vec[x]))/ncol(boot_mat)
  })
  
  model_selected <- colnames(err_mat)[which(p_vec <= alpha)]
  
  list(model_selected = model_selected, p_vec = p_vec)
}

###########

.recenter_list <- function(mat_list, vec_list){
  stopifnot(length(mat_list) == length(vec_list))
  len <- length(mat_list)
  for(i in 1:len){
    stopifnot(ncol(mat_list[[i]]) == length(vec_list[[i]]))
    
    mat_list[[i]] <- t(t(mat_list[[i]]) - vec_list[[i]])
  }
  
  mat_list
}

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

.cvc_bootstrap_trial <- function(err_mat2, mu_vec, sd_vec){
  n <- nrow(err_mat2)
  g_vec <- stats::rnorm(n)
  ((diag(g_vec) %*% err_mat2) - (g_vec %*% t(mu_vec)))/(rep(1, n) %*% t(sd_vec))
}

.cvc_bootstrap_trial_cv <- function(err_mat_list2, mu_vec_list, sd_vec){
  len <- length(err_mat_list2)
  tmp_list <- lapply(1:len, function(i){
    .cvc_bootstrap_trial(err_mat_list2[[i]], mu_vec_list[[i]], sd_vec)
  })
  
  do.call(rbind, tmp_list)
}

.extract_max <- function(vec, k){
  stopifnot(length(vec) == k * (k-1))
  sapply(1:k, function(x){
    max(vec[((x-1)*(k-1)+1) : (x*(k-1))])
  })
}

