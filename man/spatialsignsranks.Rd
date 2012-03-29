\name{Spatial signs, symmetrized signs, ranks and signed ranks}
\alias{spatial.signs}
\alias{spatial.rank}
\alias{spatial.symmsign}
\alias{spatial.signrank}
\title{Spatial signs, symmetrized signs, ranks and signed ranks}

\description{Functions to compute spatial signs, symmetrized signs, ranks and
signed ranks.}

\usage{
spatial.signs(X, center = TRUE, shape = TRUE, 
na.action = na.fail,...)
spatial.symmsign(X, shape = TRUE, na.action = na.fail, ...)
spatial.rank(X, shape = TRUE, na.action = na.fail, ...)
spatial.signrank(X, center = TRUE, shape = TRUE, 
na.action = na.fail,...)
}
\arguments{
  \item{X}{a matrix or a data frame}
  \item{center}{a vector or a logical, see details}
  \item{shape}{a matrix or a logical, see details}
  \item{\dots}{arguments that can be passed on to function used for the
 estimation of shape.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}

\details{ The spatial signs of an observed vector is simply the
vector, possibly affinely transformed first, multiplied by its
Euclidian length. See \code{\link[ICSNP]{spatial.sign}} for a precise
definition. Symmetrized spatial signs are the spatial signs of the pairwise
differences of the data 
\deqn{||x_i-x_j||^{-1}(x_i-x_j)}{||x_i-x_j||^{-1}(x_i-x_j)}
(there are \code{n} over 2 of these). Spatial
rank of an observation is the average of the signs of the differences
of that observation and the others:
\deqn{R(x_i)=ave_j\{||x_i-x_j||^{-1}(x_i-x_j)\}}{R(x_i)=ave_j\{||x_i-x_j||^{-1}(x_i-x_j)\}}
Spatial signed rank of an observation is defined as
\deqn{Q(x_i)=(R(x_i)+ave_j\{||x_i+x_j||^{-1}(x_i+x_j)\})/2}{Q(x_i)=(R(x_i)+ave_j\{||x_i+x_j||^{-1}(x_i+x_j)\})/2}
   
If a numerical value is given for \code{shape} and/or \code{center}
these are used to transform the data before the computation of signs
or ranks. A logical \code{TRUE} indicates that the shape or center should be
estimated. In this case an affine transformation that makes the
resulting signs or ranks have a covariance matrix equal or
proportional to the identity matrix and centerd on the origin is
found. A logical \code{FALSE} indicates that the null value, that is, the
identity matrix or the origin, should be used. Note that only signed
ranks depend on a center.

The value of shape and/or location used are returned as attributes.
}

\author{Seija Sirkia, \email{seija.sirkia@iki.fi}}

\references{\cite{Visuri, S., Koivunen, V. and Oja, H. (2000). Sign and rank covariance matrices. \emph{J. Statistical Planning and Inference}, 91, 557-575.}}

\seealso{\code{\link[ICSNP]{spatial.sign}} for the signs, \link[=rank.shape]{spatial sign and rank covariance matrices} and \code{\link{spatial.shape}} for the standardizing transformations}
\examples{
A<-matrix(c(1,2,-3,4),ncol=2)
X<-matrix(rnorm(100),ncol=2)\%*\%t(A)
def.par<-par(no.readonly=TRUE) # for resetting
layout(matrix(1:4,ncol=2,nrow=2,byrow=TRUE))
plot(X,col=c(2,rep(1,19)))
plot(spatial.symmsign(X),col=c(2,rep(1,19)),xlim=c(-1,1),ylim=c(-1,1))
theta<-seq(0,2*pi,length=1000)
lines(sin(theta),cos(theta))
plot(spatial.rank(X),col=c(2,rep(1,19)),xlim=c(-1,1),ylim=c(-1,1))
lines(sin(theta),cos(theta))
plot(spatial.signrank(X),col=c(2,rep(1,19)),xlim=c(-1,1),ylim=c(-1,1))
lines(sin(theta),cos(theta))
par(def.par)
}
\keyword{multivariate}
\keyword{nonparametric}
