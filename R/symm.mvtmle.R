`symm.mvtmle` <- function(X, nu=1, init=NULL, steps=Inf, eps=1e-6, maxiter=100, na.action=na.fail)
{
  X <- na.action(X)
  X <- as.matrix(X) 
  d <- dim(X)
  if(d[2]==1) return(diag(1))
 
  if (is.null(init)) init <- cov(X)
  if(is.finite(steps)) maxiter <- Inf

  iter <- 0
  V <- init
  while(TRUE){
    if(iter>=steps) return(V)
    if(iter>=maxiter) warning("maxiter reached")
    iter <- iter+1
    V.new <- matrix(.C("symm_mvtmle", as.double(X), as.double(solve(V)), 
as.integer(d), as.double(c(nu,d[2])), res=double(d[2]^2), 
PACKAGE="SpatialNP")$res, ncol=d[2], byrow=T)/(d[1]*(d[1]-1)/2)
    if(all(is.infinite(steps),mat.norm(V.new-V)<eps)){
      V <- V.new
      break
    } 
    V <- V.new 
  }
  V
}