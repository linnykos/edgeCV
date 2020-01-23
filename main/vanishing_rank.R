rm(list=ls())
library(simulation)
#library(networkCV)

set.seed(10)
ncores <- 1
doMC::registerDoMC(cores = ncores)

trials <- 50
paramMat <- as.matrix(expand.grid(round(exp(seq(log(15), log(150), length.out = 10))), c(0, 0.25, 0.5, -1),
                                  5, 200, 0.05, 0.1))
colnames(paramMat) <- c("n", "rho", "K", "trials", "alpha", "test_prop")

#############

rule <- function(vec){
  b_mat <- matrix(c(0.8, 0.7, 0.1,  0.7, 0.8, 0.1, 0.1, 0.1, 0.8), 3, 3)
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
  ecv_res <- edge_cv_sbm_sample_split(dat, k_vec = c(1:vec["K"]), test_prop = vec["test_prop"])
  err_mat <- ecv_res$err_mat
  
  cvc_res <- cvc_sbm_sample_split(err_mat, vec["trials"], vec["alpha"], verbose = F)
  
  list(err_vec = ecv_res$err_vec, p_vec = cvc_res$p_vec)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 10; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 20; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 30; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 40; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)


###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = ncores, as_list = T,
                                        filepath = "vanishing_rank_tmp.RData",
                                        verbose = T)
save.image("vanishing_rank.RData")