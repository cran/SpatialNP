\name{Shape matrices}
\alias{signrank.shape}
\alias{rank.shape}
\alias{spatial.shape}

\title{Shape matrices based on spatial ranks and signed ranks}

\description{ Iterative algorithms to find shape matrices based on
spatial signs and ranks and the k-step versions of these.  }

\usage{ 

rank.shape(X, init = NULL, steps = Inf, eps = 1e-06, maxiter = 100,
na.action = na.fail)

signrank.shape(X, location = NULL, init = NULL, steps = Inf, eps =
1e-06, maxiter = 100, na.action = na.fail)

spatial.shape(X, score = c("sign", "symmsign", "rank", "signrank"),
location = NULL, init = NULL, steps = Inf, eps = 1e-06, maxiter = 100,
na.action = na.fail)}


\arguments{
  \item{X}{a matrix or a data frame }
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{location}{an optional vector giving the location of the data} 
  \item{init}{an optional starting value for the iteration}
  \item{steps}{fixed number of iteration steps to take, if \code{Inf} iteration is repeated until convergence (or until \code{maxiter} steps)}
  \item{eps}{tolerance for convergence}
  \item{maxiter}{maximum number of iteration steps. Ignored if \code{steps} 
 is finite} 
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}  } 

\details{As \code{\link[ICSNP]{tyler.shape}} for spatial signs and
\code{\link[ICSNP]{duembgen.shape}} for spatial symmetrized signs, but for
spatial ranks and signed ranks. These are the so called inner
standardization matrices of location etc. tests based on spatial signs
and ranks. When data is standardized using these matrices the
corresponding sign or rank scores will appear \dQuote{uncorrelated}:
the corresponding outer standardization matrices will be proportional
to the identity matrix, see examples.

\code{spatial.shape} is a wrapper function for a unified access to all
four shape estimates. The choice of estimate is done via \code{score}:
\itemize{
\item \code{"sign"} for \code{\link[ICSNP]{tyler.shape}} 
\item \code{"symmsign"} for \code{\link[ICSNP]{duembgen.shape}} 
\item \code{"rank"} for \code{rank.shape} 
\item \code{"signrank"} for \code{signrank.shape}
}

\code{signrank.shape} (and \code{\link[ICSNP]{tyler.shape}}, thus also
\code{spatial.shape}) requires the location vector with respect to
which it is computed. If none is provided, vector of column means is
used.

}

\author{Seija Sirkia, \email{seija.sirkia@iki.fi}}

\seealso{\code{\link[ICSNP]{tyler.shape}}, \code{\link[ICSNP]{duembgen.shape}}, also \link[=SCov]{spatial sign and rank covariance matrices} and \link[=spatial.rank]{spatial signs and ranks}}
\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rt(150,1),ncol=3)\%*\%t(A)
signrank.shape(X)
spatial.shape(X,score="sign")
to.shape(A\%*\%t(A))
# one-step shape estimate based on spatial ranks and covariance matrix:
spatial.shape(X,score="rank",init=cov(X),steps=1)
# effect of inner standardization:
RCov(X)
RCov(X\%*\%t(chol(solve(rank.shape(X)))))
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{htest}
