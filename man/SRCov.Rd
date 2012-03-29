\name{Spatial sign and rank covariance matrices}
\alias{SCov}
\alias{SSCov}
\alias{RCov}
\alias{SRCov}
\title{Spatial sign and rank covariance matrices}
\description{
Functions to compute spatial sign, spatial symmetrized sign, spatial rank and 
spatial signed rank covariance matrices 
}
\usage{
SCov(X, location = NULL, na.action = na.fail)
SSCov(X, na.action = na.fail)
RCov(X, na.action = na.fail)
SRCov(X, location = NULL, na.action = na.fail)
}
\arguments{
  \item{X}{matrix or a data frame}
  \item{location}{numeric vector (may be missing)}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
These functions compute the matrices of the form
\deqn{
ave \{S(x_i) S^T(x_i)\}
}{ave \{S(x_i) S^T(x_i)\}}
where \eqn{S(x_i)}{S(x_i)} are the appropriate scores of the data:
spatial signs, spatial symmetrized signs, spatial ranks or spatial
signed ranks. These are the so called outer standardization matrices
of location etc. tests based on spatial signs and ranks. They are
not affine equivariant.

\code{SCov} and \code{SRCov} require a location vector with respect
to which they are computed. If none is provided, \code{SCov} uses 
spatial median and \code{SRCov} uses Hodges-Lehmann estimator.
}


\author{Seija Sirkia, \email{seija.sirkia@iki.fi}}

\references{\cite{Visuri, S., Koivunen, V. and Oja, H. (2000). Sign and rank covariance matrices. \emph{J. Statistical Planning and Inference}, 91, 557-575.}}

\seealso{\link[=spatial.rank]{spatial
signs and ranks}, \link[=rank.shape]{corresponding shape matrices}
(inner standardization matrices)} 

\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rt(150,1),ncol=3)\%*\%t(A) 
SCov(X) 
SSCov(X) 
RCov(X) 
SRCov(X)
to.shape(A\%*\%t(A),trace=1) 
} 

\keyword{multivariate}
\keyword{nonparametric}
