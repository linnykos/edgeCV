NCV.select <- function(A, max.K, cv = 3, verbose = T){
  # initialize
  dc.avg.se <- dc.avg.log <- avg.se <- avg.log <- rep(0,max.K)
  dc.avg.se[1] <- dc.avg.log[1] <- avg.se[1] <- avg.log[1] <- Inf
  dc.dev.mat <- dc.l2.mat <- sbm.dev.mat <- sbm.l2.mat <- matrix(0,cv,max.K)
  n <- nrow(A)
  sample.index <- sample.int(n)
  max.fold.num <- ceiling(n/cv)
  fold.index <- rep(1:cv,each=max.fold.num)[1:n]
  
  # start trying each value of K
  cv.index <- fold.index[sample.index]
  for(KK in 1:max.K){
    dc.l2 <- l2 <- dc.log.like <- log.like <- rep(0,cv)
    if(verbose) print(paste("Start",KK))
    
    for(k in 1:cv){
      holdout.index <- which(cv.index==k)
      train.index <- which(cv.index!=k)
      tmp.eval <- cv.evaluate(A,train.index,holdout.index,KK)
      
      log.like[k] <- tmp.eval$loglike
      
      sbm.l2.mat[k,KK] <- l2[k] <- tmp.eval$l2
      sbm.dev.mat[k,KK] <- log.like[k] <- tmp.eval$loglike
    }
    
    avg.se[KK] <- mean(l2)
    avg.log[KK] <- mean(log.like)
    
    if(verbose) print(paste("Finish ",KK,"....",sep=""))
  }
  
  # select for DCSBM and SBM, either based on minimizing average log-likelihood or average std
  dev.model <- paste("SBM",which.min(avg.log),sep="-")
  l2.model <- paste("SBM",which.min(avg.se),sep="-")
  
  list(dev = avg.log, l2 = avg.se, sbm.l2.mat = sbm.l2.mat,
       sbm.dev.mat = sbm.dev.mat, l2.model = l2.model, dev.model = dev.model)
}

cv.evaluate <- function(A, train.index, holdout.index, K, tol = 1e-6){
  n <- nrow(A)
  A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
  n.holdout <- length(holdout.index)
  n.train <- n-n.holdout
  A1 <- A.new[1:n.train,]
  A1.svd <- irlba::irlba(A1+0.001,nu=K,nv=K)
  
  if(K==1){
    A0 <- A1[1:n.train,1:n.train]
    pb <- sum(A0)/n.train^2
    if(pb < tol) pb <- tol
    if(pb > 1- tol) pb <- 1-tol
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(pb)) - sum((1-A.2[sum.index])*log(1-pb))
    l2 <- sum((A.2[sum.index]-pb)^2)
    return(list(loglike=loglike,l2=l2))
  }
  
  V <- A1.svd$v
  km <- kmeans(V,centers=K,nstart=30,iter.max=30)
  
  degrees <- colSums(A1)
  no.edge <- sum(degrees==0)
  

  B <- matrix(0,K,K)
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      N.1i <- intersect(1:n.train,which(km$cluster==i))
      N.2i <- intersect((n.train+1):n,which(km$cluster==i))
      N.1j <- intersect(1:n.train,which(km$cluster==j))
      N.2j <- intersect((n.train+1):n,which(km$cluster==j))
      B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(length(N.1i)*length(N.1j)+length(N.1j)*length(N.2i)+length(N.1i)*length(N.2j)+1)
    }
  }
  B <- B+t(B)
  
  Theta <- matrix(0,n,K)
  for(i in 1:K){
    N.1i <- intersect(1:n.train,which(km$cluster==i))
    N.2i <- intersect((n.train+1):n,which(km$cluster==i))
    B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(length(N.1i)*(length(N.1i)-1)/2+length(N.1i)*length(N.2i)+1)
    Theta[which(km$cluster==i),i] <- 1
    
  }
  
  P.hat.holdout <- Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])
  P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
  P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
  A.2 <- A.new[(n.train+1):n,(n.train+1):n]
  sum.index <- lower.tri(A.2)
  loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))
  
  l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
  return(list(loglike=loglike,l2=l2))
}
