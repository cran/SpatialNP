`SRCov` <-
function(X,location, na.action=na.fail)
{
X <- na.action(X)
if (!all(sapply(X, is.numeric))) stop("'X' must be numeric")
X<-as.matrix(X)

if(missing(location)) location<-colMeans(X)
R<-as.matrix(signranks(X))
t(R)%*%R/dim(R)[1]
}

