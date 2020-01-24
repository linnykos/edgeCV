rm(list=ls())
library(simulation)
#library(networkCV)

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
  dat <- generate_sbm(b_mat, cluster_idx, rho)
  
  dat
}

criterion <- function(dat, vec, y){
  ecv_res <- edge_cv_sbm(dat, k_vec = c(1:vec["K"]), nfold = vec["nfold"], verbose = F)
  err_mat_list <- ecv_res$err_mat_list
  
  cvc_res <- cvc_sbm(err_mat_list, vec["trials"], vec["alpha"], verbose = F)
  
  list(err_vec = ecv_res$err_vec, p_vec = cvc_res$p_vec)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 2; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)


###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = ncores, as_list = T,
                                        filepath = "sparsity_2_tmp.RData",
                                        verbose = T)
save.image("sparsity_2.RData")