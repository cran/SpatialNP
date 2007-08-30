`ae.hl.estimate` <-
function(X,init=NULL,shape=TRUE,maxiter=500,eps=1e-6,na.action=na.fail)
{
 X<-na.action(X)
 if(!all(sapply(X,is.numeric))) stop("'X' must be numeric")
 X<-as.matrix(X)

 p<-dim(X)[2]
 if(p==1) return(median(c(pairsum(X)/2,X)))

 if(is.matrix(shape))
 {
  X<-X%*%solve(mat.sqrt(shape))
  res<-as.vector(mat.sqrt(shape)%*%spatial.median(rbind(pairsum(X)/2,X)))
  attr(res,"shape")<-shape
  return(res)
 }
 else if(is.logical(shape))
 {
  if(!shape) 
  {
   res<-spatial.median(rbind(pairsum(X)/2,X))
   attr(res,"shape")<-diag(p)
   return(res)
  }
 }
 else stop("'shape' must be a matrix or logical")

 if (is.null(init))
  init <- apply(X,2,median)
 else 
  if(any(!is.vector(init),!is.numeric(init)))
   stop("'init' must be a numeric vector or NULL")
  else if (length(init)!=p) stop("'init' is of wrong dimension")
 
 X2<-rbind(X,pairsum(X)/2)
 V0<-signrank.shape(X,init)
 A0<-solve(mat.sqrt(V0))
 differ<-Inf
 iter<-0
 while(differ>eps)
 {
  if (iter==maxiter)
  {
   stop("maxiter reached without convergence")
  }
  theta.k1<-solve(A0)%*%spatial.median(X2%*%t(A0))
  V.k1<-signrank.shape(X,theta.k1)
  A.k1 <- solve(mat.sqrt(V.k1))
  theta.k2<-solve(A.k1)%*%spatial.median(X2%*%t(A.k1))
   V.k2<-signrank.shape(X,theta.k2)
  A0<-solve(mat.sqrt(V.k2))
  differ<-sqrt(sum((theta.k1-theta.k2)^2))
  iter=iter+1
 }
 res<-as.vector(theta.k2)
 attr(res,"shape")<-V.k2
 res
}

