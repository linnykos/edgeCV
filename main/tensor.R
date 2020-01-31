rm(list=ls())
library(simulation)
library(networkCV)

set.seed(10)
trials <- 150
paramMat <- as.matrix(expand.grid(100, 10, c(seq(0.3, 0.5, length.out = 5)),
                                  c(3:5), 6 ,200, 5))
colnames(paramMat) <- c("n", "p", "rho", "K", "num_model", "trials", "nfold")

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
  
  b_tensor <- array(NA, c(vec["p"], vec["K"], vec["K"]))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  dat <- networkCV::generate_tensor(b_tensor, cluster_idx, rho)
  
  dat
}

criterion <- function(dat, vec, y){
  ecv_res <- networkCV::edge_cv_sbm_tensor(dat, k_vec = c(1:vec["num_model"]), nfolds = vec["nfold"], verbose = F)
  err_mat_list <- ecv_res$err_mat_list
  
  ecv_cvc_res <- networkCV::cvc(do.call(rbind, err_mat_list), vec["trials"])
  
  list(ecv_err_vec = ecv_res$err_vec, ecv_p_vec = ecv_cvc_res)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# idx <- 3; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)


###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 25, as_list = T,
                                        filepath = "../results/tensor_tmp.RData",
                                        verbose = T)
save.image("../results/tensor.RData")