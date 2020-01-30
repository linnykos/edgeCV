rm(list=ls())
source("../experiment/ncv_source.R")
library(simulation)
library(networkCV)


set.seed(10)
trials <- 100
paramMat <- as.matrix(expand.grid(c(600), c(0.05, 0.1, 0.2),
                                  c(3:5),
                                  6 , 200, 5))
colnames(paramMat) <- c("n", "rho", "K", "num_model", "trials", "nfold")

#############

rule <- function(vec){
  b_mat <- matrix(1, vec["K"], vec["K"]) + 2*diag(vec["K"])
  n <- vec["n"]
  n_each <- round(n/vec["K"])
  cluster_idx <- rep(1:(vec["K"]-1), each = n_each)
  cluster_idx <- c(cluster_idx, rep(vec["K"], n - length(cluster_idx)))
  dat <- networkCV::generate_sbm(b_mat, cluster_idx, vec["rho"])
  
  dat
}

criterion <- function(dat, vec, y){
  cat(paste0(".", y, "."))
  ncv_res <- NCV(dat, K.vec = c(1:vec["num_model"]), n.fold = vec["nfold"])
  err_mat<- ncv_res$err_mat
  
  cvc_res <- networkCV::cvc(err_mat, vec["trials"])
  
  list(err_vec = ncv_res$err_vec, p_vec = cvc_res)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)

########


res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 20, as_list = T,
                                        filepath = "../results/ncv_tmp.RData",
                                        verbose = T)
save.image("../results/ncv.RData")
