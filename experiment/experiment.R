rm(list=ls())
library(networkCV)
set.seed(10)
b_mat_truth <- matrix(c(1/4, 1/2, 1/4,  1/2, 1/4, 1/4,  1/4, 1/4, 1/6), 3, 3)
cluster_idx_truth <- c(rep(1, 10), rep(2, 10), rep(3, 10))
dat <- networkCV::generate_sbm(b_mat_truth, cluster_idx_truth)
k <- 5
err_mat_list <- networkCV::edge_cv_sbm(dat, k_vec = c(1:k), nfold = 5, verbose = F)$err_mat_list
trials <- 10
alpha <- 0.05
verbose = F
ncores = 5

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
  networkCV:::.clean_err_mat(err_mat)
})

mu_vec_list <- lapply(err_mat_list2, function(err_mat2){colMeans(err_mat2)})

err_mat_list2_recentered <- networkCV:::.recenter_list(err_mat_list2, mu_vec_list)

mu_vec <- colMeans(do.call(rbind, mu_vec_list))
sd_vec <- apply(do.call(rbind, err_mat_list2_recentered), 2, stats::sd)

# compute test statistics
test_vec <- networkCV:::.extract_max(sqrt(n)*mu_vec/sd_vec, k)

# compute null distribution
func <- function(b){
  if(verbose && b %% floor(trials/10) == 0) cat('*')
  
  set.seed(sum(sd_vec)*100+10*b)
  tmp_mat <- networkCV:::.cvc_bootstrap_trial_cv(err_mat_list2, mu_vec_list, sd_vec)
  networkCV:::.extract_max((1/sqrt(n)) * colSums(tmp_mat), k)
}

b <- 0 #debugging purposes
boot_mat <- foreach::"%dopar%"(foreach::foreach(b = 1:trials), func(b))
boot_mat <- do.call(rbind, boot_mat)

boot_mat2 <- sapply(1:trials, func)
