\name{Spatial location}
\alias{ae.hl.estimate}
\alias{spatial.location}
\title{Multivariate location estimates based on spatial signs and signed ranks}

\description{ Spatial median, multivariate Hodges-Lehmann estimate of
location and their affine equivariant versions }

\usage{

ae.hl.estimate(X, init=NULL, shape=TRUE, maxiter = 500, eps = 1e-6,
na.action = na.fail)

spatial.location(X, score = c("sign", "signrank"), init = NULL, shape
= TRUE, maxiter = 500, eps = 1e-6, na.action = na.fail) }

\arguments{
  \item{X}{a matrix or a data frame}
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{init}{an optional vector giving the initial point of the iteration}
  \item{shape}{logical, or a matrix. See details}
  \item{eps}{tolerance for convergence}
  \item{maxiter}{maximum number of iteration steps}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
 } 

\details{Hodges-Lehmann estimate of multivariate location is the
spatial median (see \code{\link[ICSNP]{spatial.median}}) of the observations
and their pairwise averages. Since spatial median is not affine
equivariant, neither is multivariate Hodges-Lehmann estimate. Affine
equivariance can be achieved by simultaneously estimating the
corresponding shape, as proposed for the spatial median by
Hettmansperger and Randles (2002). For spatial median the corresponding
shape is \code{\link[ICSNP]{tyler.shape}} and the simultaneous estimation is
implemented as \code{\link[ICSNP]{HR.Mest}}, see its documentation for
further details on the algorithm. For the Hodges-Lehmann estimate it
is \code{\link{rank.shape}}.

\code{spatial.location} is a wrapper function for a unified access to
both location estimates. The choice of estimate is done via
\code{score}:
\item \code{"sign"} for spatial median 
\item \code{"signrank"} for Hodges-Lehmann estimate 

If a matrix (must be symmetric and positive definite, but this is not
checked) is given as \code{shape} the location estimate is found with
respect to that shape and no further shape estimation is done. If a
logical \code{TRUE} is given as \code{shape} the shape is estimated
and consequently the affine equivariant version of the location
estimate is found.  If \code{shape} is \code{FALSE} then shape
estimation is not done and the non affine equivariant version of the
location estimate is found.

}  

\value{
The estimate vector with the (final estimate of or given) shape matrix 
as attribute \code{"shape"}.
}
\references{Hettmansperger, T. and Randles, R. (2002) A Practical Affine Equivariant Multivariate Median, Biometrika, 89, pp. 851-860}
\author{Seija Sirkia, \email{ssirkia@maths.jyu.fi}}

\seealso{\code{\link[ICSNP]{spatial.median}}, \code{\link[ICSNP]{HR.Mest}} \code{\link[ICSNP]{signrank.shape}}}

\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rt(150,1),ncol=3)\%*\%t(A)
spatial.location(X,score="signrank")
spatial.location(X,score="sign")
#compare with:
colMeans(X)
ae.hl.estimate(X,shape=A\%*\%t(A))
ae.hl.estimate(X,shape=FALSE)
}
\keyword{multivariate}
\keyword{nonparametric}
