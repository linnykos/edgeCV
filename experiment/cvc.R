CleanErrMat <- function(err.mat) {
  # This is a subroutine called in the main function "CVPerm"
  # It removes identical columns in the cross-validated loss matrix
  n2 <- nrow(err.mat)
  M <- ncol(err.mat)
  err.mean <- apply(err.mat, 2, mean)
  err.mat.center <- err.mat - matrix(err.mean, nrow = n2, ncol = M, byrow = T)
  
  res <- list()
  list.length <- 0
  remaining.index <- 1:M
  
  while (length(remaining.index) > 0) {
    m <- remaining.index[1]
    err.diff.center <- err.mat.center[,m] - err.mat.center
    err.mean.diff <- err.mean[m] - err.mean
    sd.vec <- apply(err.diff.center, 2, sd)
    J.1 <- which(sd.vec==0 & err.mean.diff > 0)
    J.2 <- which(sd.vec==0 & err.mean.diff == 0)
    J.3 <- which(sd.vec==0 & err.mean.diff < 0)
    if (length(J.1) == 0) {
      list.length <- list.length+1
      res[[list.length]] <- J.2
    }
    remaining.index <- setdiff(remaining.index, union(J.2,J.3))
  }
  return(list(err.mat.c = err.mat[,sapply(res,function(x){x[1]})], ind=res))
}

CVPerm <- function(err.mat, B = 200, screen=T, alpha.s=0.005) {
  # This is the main function of cross-validation with confidence.
  # Requires subroutine "CleanErrMat"
  # Input: "err.mat" is a matrix where each row corresponds to a data point
  #            and each column corresponds to a candidate model.
  #            err.mat[i,j] records the cross-validated loss evaluated at
  #            the i'th data point and the j'th candidate model.
  #            For example, in least square linear regression this is the
  #            cross-validated squared residual 
  #        "B" is the number of bootstrap samples
  #        "screen" is an indicator if pre-screening is used for the test
  #        "alpha.s" is the threshold used in pre-screening
  clean.res <- CleanErrMat(err.mat)
  err.mat.c <- clean.res$err.mat.c
  err.mat.c.ind <- clean.res$ind
  n2 <- nrow(err.mat.c)
  M <- ncol(err.mat.c)
  err.mean <- apply(err.mat.c, 2, mean)
  err.mat.center <- err.mat.c - matrix(err.mean, nrow = n2, ncol = M, byrow = T)
  sign.mat <- matrix(rnorm(n2*B),ncol = B)
  norm.quantile <- qnorm(1-alpha.s/(M-1))
  screen.th <- ifelse(norm.quantile^2>=n2,Inf,norm.quantile / sqrt(1-norm.quantile^2/n2))
  sgmb.p.val <- sapply(1:M, function(m){
    err.diff.center <- err.mat.center[,m] - err.mat.center[, setdiff(1:M, m),drop=F]
    err.mean.diff <- err.mean[m] - err.mean[setdiff(1:M, m)]
    sd.vec <- apply(err.diff.center, 2, sd)
    err.mean.diff.scale <- err.mean.diff / sd.vec
    test.stat <- sqrt(n2)*max(err.mean.diff.scale)
    if (test.stat >= screen.th) {
      return(alpha.s)
    }
    if (screen) {        
      J.screen <- which(sqrt(n2) * err.mean.diff.scale > -2*screen.th)
      if (length(J.screen)==0) {return(1-alpha.s)}
    } else {
      J.screen <- 1:(M-1)
    }
    err.diff.center.scale <- err.diff.center[,J.screen] / 
      matrix(sd.vec[J.screen], nrow = n2, ncol = length(J.screen), byrow=T)
    test.stat.vec <- sapply(1:B,
                            function(ib){max(apply(err.diff.center.scale*sign.mat[,ib],2,mean))})
    return(mean(test.stat.vec>max(err.mean.diff.scale[J.screen])))
  })
  res <- rep(0,ncol(err.mat))
  for (i in 1:length(err.mat.c.ind)) {
    res[err.mat.c.ind[[i]]] <- sgmb.p.val[i]
  }
  return(res)
}