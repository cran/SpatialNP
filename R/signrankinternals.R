norm<-function(X) 
.C("norming", as.double(X), as.integer(dim(X)), res=double(dim(X)[1]),PACKAGE="ICSNP")$res

pairdiff<-function(X)
{
d<-dim(X)
matrix(.C("pairdiff", as.double(X),as.integer(d), res=double(choose(d[1],2)*d[2]),PACKAGE="ICSNP")$res,ncol=d[2],byrow=T)
}

pairsum<-function(X)
{
d<-dim(X)
matrix(.C("pairsum", as.double(X),as.integer(d), res=double(choose(d[1],2)*d[2]),PACKAGE="ICSNP")$res,ncol=d[2],byrow=T)
}

sumsignout<-function(X)
{
#have to be careful with zero divisions:
ind<-numeric(0)
ind<-apply(X,1,setequal,y=0)
if(length(ind)>0) 
 X<-X[!ind,]

d<-dim(X)
tmp<-matrix(.C("sum_of_sign_outers", as.double(X),as.integer(d), res=double(d[2]^2),PACKAGE="ICSNP")$res,ncol=d[2],byrow=T)
}

ranks<-function(X)
{
d<-dim(X)
if(any(d[2]==1,is.vector(X))) return(rank(X))
tmp<-matrix(.C("spatial_ranks", as.double(X),as.integer(d), res=double(d[1]*d[2]),PACKAGE="ICSNP")$res,ncol=d[2],byrow=T)

#need to check whether multiple values caused a problem:
ind<-which(!is.finite(tmp[,1]))
if(length(ind)==0) return(tmp)
#else the problem ranks need to be recalculated:
for(i in ind) {
 s<-sweep(X,2,X[i,]) 
 r<-norm(s)
 r[r==0]<-1
 tmp[i,]<-apply(sweep(s,1,r,"/"),2,mean)
}
tmp
}

signranks<-function(X)
{
d<-dim(X)
tmp<-matrix(.C("signed_ranks", as.double(X),as.integer(d), res=double(d[1]*d[2]),PACKAGE="ICSNP")$res,ncol=d[2],byrow=T)
#as in ranks, possibly need to recalculate stuff
ind<-which(!is.finite(tmp[,1]))
if(length(ind)==0) return(tmp)
#else 
for(i in ind) {
 sm<-sweep(X,2,X[i,]) 
 sp<-sweep(X,2,X[i,],"+")
 rm<-norm(sm)
 rp<-norm(sp)
 rm[rm==0]<-1
 rp[rp==0]<-1
 tmp[i,]<-(apply(sweep(sm,1,rm,"/"),2,mean)+apply(sweep(sp,1,rp,"/"),2,mean))/2
}
tmp
}

Q2internal<-function(X)
{
d<-dim(X)
.C("Q2internals",as.double(X),as.integer(d), res=double(d[2]^2+d[2]^4),PACKAGE="ICSNP")$res
}

gen.inv<-function(M)
{
p<-sqrt(dim(M)[1])
eig<-eigen(M)
eig$values[1:((p+2)*(p-1)/2)]<-eig$values[1:((p+2)*(p-1)/2)]^-1
res<-eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
Re(res) #because it really is real!
}

Cpp<-function(p)
{
I<-diag(p)
vecI<-I
dim(vecI)<-NULL
J<-vecI%*%t(vecI)
K<-matrix(0,ncol=p^2,nrow=p^2)
for (i in 1:p) for (j in 1:p) K<-K+kronecker(outer(I[,i],I[,j]),outer(I[,j],I[,i]))
(diag(p^2)+K)/2-J/p
}

covshape<-function(x) to.shape(cov(x))

mat.sqrt<-function(A)
# Returns the square root matrix of the given matrix.
{
 eig<-eigen(A)
 eig$vectors%*%(diag(eig$values))^(1/2)%*%t(eig$vectors)
}

mat.norm<-function(A)
# Returns the matrix norm of the given matrix.
{
 sqrt(sum(diag(t(A)%*%A)))
}
