rm(list=ls())
library(simulation)
library(networkCV)

set.seed(10)
trials <- 100
paramMat <- as.matrix(expand.grid(c(600), c(seq(0.1, 0.6, length.out = 13)),
                                  c(3:5),
                                  6 , 200, 5))
colnames(paramMat) <- c("n", "rho", "K", "num_model", "trials", "nfold")

#############

rule <- function(vec){
  b_mat <- matrix(0.2, vec["K"], vec["K"]) + 0.6*diag(vec["K"])
  n <- vec["n"]
  n_each <- round(n/vec["K"])
  cluster_idx <- rep(1:(vec["K"]-1), each = n_each)
  cluster_idx <- c(cluster_idx, rep(vec["K"], n - length(cluster_idx)))
  if(vec["rho"] >= 0){
    rho <- 1/(n^vec["rho"])
  } else {
    rho <- log(n)/n
  }
  dat <- networkCV::generate_sbm(b_mat, cluster_idx, rho)
  
  dat
}

criterion <- function(dat, vec, y){
  ecv_res <- networkCV::edge_cv_sbm(dat, k_vec = c(1:vec["num_model"]), nfold = vec["nfold"], verbose = F)
  err_mat_list <- ecv_res$err_mat_list
  
  cvc_res <- networkCV::cvc(do.call(rbind, err_mat_list), vec["trials"])
  
  list(err_vec = ecv_res$err_vec, p_vec = cvc_res)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 24; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)


###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 20, as_list = T,
                                        filepath = "../results/sparsity_4_tmp.RData",
                                        verbose = T)
save.image("../results/sparsity_4.RData")