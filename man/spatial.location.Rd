\name{Spatial location}
\alias{ae.hl.estimate}
\alias{ae.spatial.median}
\alias{spatial.location}
\title{Multivariate location estimates based on spatial signs and signed ranks}

\description{ Iterative algorithms to find spatial median, multivariate Hodges-Lehmann estimate of
location, their affine equivariant versions and k-step versions of these. }

\usage{

spatial.location(X, score = c("sign", "signrank"), init = NULL, 
shape = TRUE, steps = Inf, maxiter = 500, eps = 1e-6, 
na.action = na.fail) 

ae.spatial.median(X, init = NULL, shape = TRUE, steps = Inf, 
maxiter = 500, eps = 1e-6, na.action = na.fail)

ae.hl.estimate(X, init = NULL, shape = TRUE, steps = Inf, 
maxiter = 500, eps = 1e-06, na.action = na.fail)
}

\arguments{
  \item{X}{a matrix or a data frame}
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{init}{an optional vector giving the initial point of the iteration}
  \item{shape}{logical, or a matrix. See details}
  \item{steps}{fixed number of iteration steps to take, if \code{Inf} iteration is repeated until convergence (or until \code{maxiter} steps)}
  \item{eps}{tolerance for convergence}
  \item{maxiter}{maximum number of iteration steps}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
 } 

\details{Spatial median and Hodges-Lehmann estimator (spatial median of the pairwise differences) are not affine equivariant. Affine
equivariance can be achieved by simultaneously estimating the
corresponding shape, as proposed for the spatial median by
Hettmansperger and Randles (2002). For spatial median the corresponding
shape is \code{\link{signs.shape}} and for the Hodges-Lehmann estimate it
is \code{\link{signrank.shape}}.

\code{spatial.location} is a wrapper function for a unified access to
both location estimates. The choice of estimate is done via
\code{score}:
\itemize{
\item \code{"sign"} for spatial median 
\item \code{"signrank"} for Hodges-Lehmann estimate 
}

If a matrix (must be symmetric and positive definite, but this is not
checked) is given as \code{shape} the location estimate is found with
respect to that shape and no further shape estimation is done. If a
logical \code{TRUE} is given as \code{shape} the shape is estimated
and consequently the affine equivariant version of the location
estimate is found.  If \code{shape} is \code{FALSE} then shape
estimation is not done and the non affine equivariant versions of the
location estimate, that is the spatial median and the Hodges-Lehmann estimate are found.

}  

\value{
The estimate vector with the (final estimate of or given) shape matrix 
as attribute \code{"shape"}.
}
\references{Hettmansperger, T. and Randles, R. (2002) A Practical Affine Equivariant Multivariate Median, Biometrika, 89, pp. 851-860}

\author{Seija Sirkia, \email{seija.sirkia@iki.fi}, Jari Miettinen, \email{jari.p.miettinen@aalto.fi}}

\seealso{\code{\link[ICSNP]{spatial.median}}, \code{\link{signrank.shape}}}

\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rnorm(3000),ncol=3)\%*\%t(A)
spatial.location(X,score="signrank")
spatial.location(X,score="sign")
#compare with:
colMeans(X)
ae.hl.estimate(X,shape=A\%*\%t(A))
ae.hl.estimate(X,shape=FALSE)
}
\keyword{multivariate}
\keyword{nonparametric}
