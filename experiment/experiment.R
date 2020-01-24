rm(list=ls())
library(simulation)
library(networkCV)

set.seed(10)
ncores <- NA
doMC::registerDoMC(cores = ncores)

trials <- 50
paramMat <- as.matrix(expand.grid(c(30, 100, 300, 1000), c(0, 0.25, 0.5),
                                  5, 200, 0.05, 5))
colnames(paramMat) <- c("n", "rho", "K", "trials", "alpha", "nfold")

#############

rule <- function(vec){
  b_mat <- matrix(0.2, 3, 3) + 0.6*diag(3)
  n <- vec["n"]
  n3 <- round(n/3)
  cluster_idx <- c(rep(1, n3), rep(2, n3), rep(3, vec["n"]-2*n3))
  if(vec["rho"] >= 0){
    rho <- 1/(n^vec["rho"])
  } else {
    rho <- log(n)/n
  }
  dat <- networkCV::generate_sbm(b_mat, cluster_idx, rho)
  
  dat
}

########

idx <- 3; y <- 1; set.seed(y)
vec <- paramMat[idx,]
dat <- rule(vec)

ecv_res <- networkCV::edge_cv_sbm(dat, k_vec = c(1:vec["K"]), nfold = vec["nfold"], verbose = F)
err_mat_list <- ecv_res$err_mat_list

# cvc_res <- networkCV::cvc_sbm(err_mat_list, vec["trials"], vec["alpha"], verbose = F, ncores = 20)
trials = vec["trials"]
alpha = vec["alpha"]
verbose = T
ncores = 20

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
  tmp <- networkCV:::.extract_max((1/sqrt(n)) * colSums(tmp_mat), k)
  
  if(verbose) print(paste0("Trial ", b, ": Values ", paste0(round(tmp,2), collapse = ", ")))
  tmp
}

b <- 0 #debugging purposes
boot_mat <- foreach::"%dopar%"(foreach::foreach(b = 1:trials), func(b))
