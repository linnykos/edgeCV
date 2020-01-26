args <- 7

library(networkCV)

set.seed(10)
trials <- 3
paramMat <- as.matrix(expand.grid(c(600), c(0.5),
                                  5, 200, 0.05, 5))
colnames(paramMat) <- c("n", "rho", "K", "trials", "alpha", "nfold")

##########

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
  
  cvc_res <- networkCV::cvc_sbm(err_mat_list, vec["trials"], vec["alpha"], verbose = T, ncores = NA)
  
  print(ecv_res$err_vec)
  print(cvc_res$p_vec)
  
  list(err_vec = ecv_res$err_vec, p_vec = cvc_res$p_vec)
}

##############

idx <- 1
res_list <- vector("list", trials)
save.image(paste0("sparsity_largescale_", args, ".RData"))

for(i in 1:trials){
  y <- args*20+i
  set.seed(y)
  
  res_list[[i]] <- criterion(rule(paramMat[idx,]), paramMat[idx,], y)
  save.image(paste0("sparsity_largescale_", args, ".RData"))
}

save.image(paste0("sparsity_largescale_", args, ".RData"))
