`RCov` <-
function(X,na.action=na.fail)
{
X <- na.action(X)
if (!all(sapply(X, is.numeric))) stop("'X' must be numeric")
X<-as.matrix(X)

R<-as.matrix(ranks(X))
t(R)%*%R/dim(R)[1]
}

