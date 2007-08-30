sr.regression<- function(formula,data=NULL,score=c("sign","rank"),ae=TRUE,eps=1E-6,na.action=na.fail)
{
  score<-match.arg(score)

  mydata<-srreg.formula(formula,data)
  mydata$X<-na.action(mydata$X)
  if(is.vector(mydata$Y)) mydata$Y<-as.matrix(mydata$Y)
  if(!is.null(foo<-attr(mydata$X,"na.action"))) mydata$Y<-mydata$Y[-foo,]
  mydata$Y<-na.action(mydata$Y)
  if(!is.null(foo<-attr(mydata$Y,"na.action"))) mydata$X<-mydata$X[-foo,]


  if(score=="rank")
# Recursive! Regression based on ranks is in fact
# regression based on signs of the pairwise differences
  {
   int=FALSE
   if( length(i<-which(attr(mydata$X,"assign")==0)) > 0 ) {
    int=TRUE
   if(length(i) == dim(mydata$X)[2]) {
    temp<-rbind(ae.hl.estimate(mydata$Y,shape=TRUE,eps=eps))
    rownames(temp)<-c("(Intercept)")
    return(temp)
   }
   Xd<-pairdiff(mydata$X[,-i])
   }
   else Xd<-pairdiff(mydata$X)
   Yd <- pairdiff(mydata$Y)
  
   difdata <- data.frame(Y=Yd,X=Xd)
  
   B <- sr.regression(Yd~Xd-1,difdata,score="sign",ae,eps)

   # possibly compute the intercept #
   if(int) {
   Ehat <- mydata$Y - mydata$X[,-1] %*% B
   b0 <- ae.hl.estimate(Ehat,shape=TRUE,eps=eps)
   }
   else b0 <- rep(NA,dim(mydata$Y)[2])
  
   temp<-rbind(b0,B)
   if(int) rownames(temp)<-colnames(mydata$X)
   else rownames(temp)<-c("(Intercept)", colnames(mydata$X))
   return(temp)
  } #end rank based
  
  X <- matrix(mydata$X,nrow=nrow(mydata$X))
  Y <- matrix(mydata$Y,nrow=nrow(mydata$Y))
  
  n<-dim(X)[1]
  p<-dim(X)[2]
  d<-dim(Y)[2]

  if(d==1) ae=FALSE

  # Preliminary estimate/starting value #
  B.new <- solve(t(X)%*%X)%*%t(X)%*%Y

  if (ae==TRUE)
  {
    # estimate residuals by a preliminary estimate #
    Ehat <- Y - X %*% B.new
    V <- tyler.shape(Ehat)
    sqrt.V <- mat.sqrt(V)
    Y <- t(solve(sqrt.V) %*% t(Y))
    B.new <- solve(t(X)%*%X)%*%t(X)%*%Y
  }
  
  B <- matrix(c(rep(0,p*d)),nrow=p,ncol=d)
    
  while (sum((B-B.new)**2)>eps)
  {
    B <- B.new
    
    # estimated residuals #
    Ehat <- Y - X %*% B
  
    # Euclidian norms of the residuals #
    r <- norm(Ehat)
    
    # investigate if there are at least p+1 zero residuals, and if so, #
    # investigate the value of the objective function to see if this   #
    # is the solution #
    if (max(sort(r)[1:(p+1)])<eps)
    {
      index<-(1:n)[r>eps]
      if (sum(apply(Ehat[index,]/r[index],2,sum)**2)>min(1E-6,eps)) 
        B.new <- solve(t(X/(cbind(r)%*%rep(1,p)))%*%X)%*%(t(X)%*%(Y/(cbind(r)%*%rep(1,d))))
    } 
    else
      B.new <- solve(t(X/(cbind(r)%*%rep(1,p)))%*%X)%*%(t(X)%*%(Y/(cbind(r)%*%rep(1,d))))
  }

  # Transform back if one-step T-R procedure is used #
   if (ae==TRUE) B.new <- t(sqrt.V %*% t(B.new))

  rownames(B.new)<-colnames(mydata$X)
  B.new
}
