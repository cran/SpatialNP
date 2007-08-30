`SCov` <-
function(X,location,na.action=na.fail)
{
X <- na.action(X)
if (!all(sapply(X, is.numeric))) stop("'X' must be numeric")
X<-as.matrix(X)

if(missing(location)) location<-colMeans(X)
n<-dim(X)[1]
if (dim(X)[2]<2) return(diag(1))
sumsignout(sweep(X,2,location))/n
}

