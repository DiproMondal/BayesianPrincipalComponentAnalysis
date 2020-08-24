set.seed(4)
X <- matrix(sample(seq(1,50),size=100, replace=T), ncol=5, nrow=20, byrow=T)
#write.table(X, file="mymatrix.txt", row.names=FALSE, col.names=FALSE)

ppca <- function(X,d = 2, scale =TRUE){
  if (d >= (dim(X)[2])){
    stop("Check number of reduced dimension")
  }
  tryCatch({
  if (scale == TRUE){
    X <- (apply(X, 2, function(x) (x-min(x))))#/(max(x)-min(x))))
  }
  center_X <- t(apply(X, 1, function(x) x-rowMeans(t(X))))
  covar <- cov(center_X)
  eigs<- eigen(covar)
  eigvals_d <- eigs$values[1:d]
  eigvecs_d <- eigs$vectors[,1:d]
  
  Lambda_d <- diag(x=eigvals_d)
  sigmahat <- sum(eigs$values[-c(1:d)])/(dim(X)[2]-d)
  Id <- diag(nrow=d)
  W <- eigvecs_d %*% (sqrt(Lambda_d - sigmahat^2 * Id))
  Mhat <- t(W) %*% (W) + sigmahat^2 * Id
  Means <- center_X %*% W %*% solve(Mhat)
  P <- t(eigvecs_d) %*% t(center_X)
  ppca_out<- t(P)
  
  return(list('principal_components'=ppca_out, 'posterior_means'=Means, 
              'residual_variance'=sigmahat))
  },error= function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ppca(X,3)