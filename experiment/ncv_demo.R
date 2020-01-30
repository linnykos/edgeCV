rm(list=ls())
source("../experiment/ncv_source.R")

vec <- c(n = 100, rho = 0.2, K = 4, num_model = 6, trials = 200, nfold = 5)

# generate data
set.seed(98)
b_mat <- matrix(0.2, vec["K"], vec["K"]) + 0.6*diag(vec["K"])
n <- vec["n"]
n_each <- round(n/vec["K"])
cluster_idx <- rep(1:(vec["K"]-1), each = n_each)
cluster_idx <- c(cluster_idx, rep(vec["K"], n - length(cluster_idx)))
rho <- 1/(n^vec["rho"])
dat <- generate_sbm(b_mat, cluster_idx, rho)

############

ncv_res <- NCV(dat)