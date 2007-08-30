`SSCov` <-
function(X,na.action=na.fail)
{
X <- na.action(X)
if (!all(sapply(X, is.numeric))) stop("'X' must be numeric")
X<-as.matrix(X)
if (dim(X)[2]<2) return(diag(1))

tmp<-pairdiff(X)
tmp2<-sumsignout(tmp)/dim(tmp)[1]
tr<-sum(diag(tmp2))
if (tr!=1 & tr!=0)
tmp2<-to.shape(tmp2,trace=1)
tmp2
}

