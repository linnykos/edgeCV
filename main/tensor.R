rm(list=ls())
library(simulation)
library(networkCV)

set.seed(10)
trials <- 50
paramMat <- as.matrix(expand.grid(c(10, 30, 100, 200), c(0, 0.25, 0.5),
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
  
  b_tensor <- array(NA, c(vec["n"], 3, 3))
  for(i in 1:dim(b_tensor)[1]){
    b_tensor[i,,] <- b_mat
  }
  
  dat <- networkCV::generate_tensor(b_tensor, cluster_idx, rho)
  
  dat
}

criterion <- function(dat, vec, y){
  print(y)
  ecv_res <- networkCV::edge_cv_sbm_tensor(dat, k_vec = c(1:vec["K"]), trials = 5, test_prop = 0.1, verbose = F)
  
  list(err_vec = ecv_res$err_vec)
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