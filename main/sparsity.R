rm(list=ls())
library(simulation)
library(networkCV)

set.seed(10)
trials <- 100
paramMat <- as.matrix(expand.grid(c(100), c(seq(0, 0.5, length.out = 11), -1),
                                  5, 200, 5))
colnames(paramMat) <- c("n", "rho", "K", "trials", "nfold")

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

criterion <- function(dat, vec, y){
  print(y)
  ecv_res <- networkCV::edge_cv_sbm(dat, k_vec = c(1:vec["K"]), nfold = vec["nfold"], verbose = F)
  err_mat_list <- ecv_res$err_mat_list
  
  cvc_res <- networkCV::cvc(do.call(rbind, err_mat_list), vec["trials"])
  
  list(err_vec = ecv_res$err_vec, p_vec = cvc_res)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 5; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)


###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 20, as_list = T,
                                        filepath = "sparsity_3_tmp.RData",
                                        verbose = T)
save.image("sparsity_3.RData")