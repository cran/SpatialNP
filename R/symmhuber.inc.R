`symmhuber.inc` <- function(X, qg=0.9, m=10, init=NULL, steps=Inf, permute=TRUE, eps=1e-6, maxiter=100, na.action = na.fail)
{
  X <- na.action(X)
  X <- as.matrix(X) 
  d <- dim(X) 
  if(d[2]==1) return(diag(1))
 
  if (is.null(init)) init <- cov(X)
  if(is.finite(steps)) maxiter <- Inf

  c.square <- 2*qchisq(qg, d[2])
  sigma.square <- 2*pchisq(c.square/2, d[2]+2)+(c.square/d[2])*(1 - qg)

  if(permute) X <- X[sample(1:d[1]),]
  
  iter <- 0
  V <- init
  while(TRUE){
    if(iter>=steps) return(V)
    if(iter>=maxiter) warning("maxiter reached")
    iter <- iter+1
    sqrtV <- mat.sqrt(V)
    V.new <- matrix(.C("symm_huber_inc", as.double(X), as.double(solve(V)), 
as.integer(d), as.double(c.square), as.double(sigma.square), as.integer(m), 
res=double(d[2]^2), PACKAGE="SpatialNP")$res, ncol=d[2], 
byrow=T)/(m*(2*d[1]-m-1)/2)

    if(all(is.infinite(steps),mat.norm(V.new-V)<eps)){
      V <- V.new
      break
    } 
    V <- V.new 
  }
  V
}


