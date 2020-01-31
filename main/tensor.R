rm(list=ls())
library(simulation)
library(networkCV)

set.seed(10)
trials <- 20
paramMat <- as.matrix(expand.grid(100, 10, c(seq(0, 0.5, length.out = 6)),
                                  5, 200, 5))
colnames(paramMat) <- c("n", "p", "rho", "K", "trials", "nfold")

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
  
  b_tensor <- array(NA, c(vec["p"], 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  dat <- networkCV::generate_tensor(b_tensor, cluster_idx, rho)
  
  dat
}

criterion <- function(dat, vec, y){
  ecv_res <- networkCV::edge_cv_sbm_tensor(dat, k_vec = c(1:vec["K"]), nfolds = vec["nfold"], verbose = F)
  err_mat_list <- ecv_res$err_mat_list
  
  ecv_cvc_res <- networkCV::cvc(do.call(rbind, err_mat_list), vec["trials"])
  
  list(ecv_err_vec = ecv_res$err_vec, ecv_p_vec = ecv_cvc_res)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 3; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)


###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 20, as_list = T,
                                        filepath = "tensor_tmp.RData",
                                        verbose = T)
save.image("tensor.RData")