rm(list=ls())

# source the relevant files
source("../networkCV/R/generate_sbm.R")
source("../networkCV/R/edge_cv.R")
source("../networkCV/R/cvc.R")

# set up the parameters needed for this demo
# n: number of nodes (total)
# rho_n: sparsity level, where 0.2 means rho_n = 1/n^(0.2)
# K: number of true clusters
# num_model: number of models to test (i.e., 6 means we try 1 through 6 number of clusters)
# trials: number of bootstrap trials per p-value
# nfold: number of folds
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

# run ECV and CVC
ecv_res <- edge_cv_sbm(dat, k_vec = c(1:vec["num_model"]), nfolds = vec["nfold"], verbose = F)
err_mat_list <- ecv_res$err_mat_list

cvc_res <- cvc(do.call(rbind, err_mat_list), vec["trials"])

# look at ecv results
ecv_res$err_vec #vector of mean test errors, one per model
which.min(ecv_res$err_vec) == vec["K"]

cvc_res #vector of p-values, one per model
# the above results are incoherent, as models with Khat = 2,5,6 were selected...
# i'm not sure if this implies the code is incorrect??

#############################

# let's do a simulation where things work out nicely (an easy setting)
vec <- c(n = 100, rho = 0, K = 3, num_model = 6, trials = 200, nfold = 5)

# generate data
set.seed(1)
b_mat <- matrix(0.2, vec["K"], vec["K"]) + 0.6*diag(vec["K"])
n <- vec["n"]
n_each <- round(n/vec["K"])
cluster_idx <- rep(1:(vec["K"]-1), each = n_each)
cluster_idx <- c(cluster_idx, rep(vec["K"], n - length(cluster_idx)))
rho <- 1/(n^vec["rho"])
dat <- generate_sbm(b_mat, cluster_idx, rho)

# run ECV and CVC
ecv_res <- edge_cv_sbm(dat, k_vec = c(1:vec["num_model"]), nfolds = vec["nfold"], verbose = F)
err_mat_list <- ecv_res$err_mat_list

cvc_res <- cvc(do.call(rbind, err_mat_list), vec["trials"])

which.min(ecv_res$err_vec) == vec["K"]
min(which(cvc_res > 0.05)) == vec["K"]

