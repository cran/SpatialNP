\name{Shape matrices}
\alias{rank.shape}
\alias{signs.shape}
\alias{signrank.shape}
\alias{spatial.shape}
\alias{symmsign.shape}
\alias{symmhuber.shape}

\title{Shape matrices based on spatial ranks and signed ranks}

\description{ Iterative algorithms to find shape matrices based on
spatial signs and ranks and the k-step versions of these.  }

\usage{ 

spatial.shape(X, score = c("sign", "symmsign", "rank", "signrank"), 
fixed.loc = FALSE, location = NULL, init = NULL, steps = Inf, 
eps = 1e-06, maxiter = 100, na.action = na.fail)

signs.shape(X, fixed.loc = FALSE, location = NULL, init = NULL, 
steps = Inf, eps = 1e-6, maxiter = 100, na.action = na.fail) 

symmsign.shape(X, init = NULL, steps = Inf, eps = 1e-6, 
maxiter = 100, na.action = na.fail)

rank.shape(X, init = NULL, steps = Inf, eps = 1e-06, 
maxiter = 100, na.action = na.fail)

signrank.shape(X, fixed.loc = FALSE, location = NULL, init = NULL, 
steps = Inf, eps = 1e-06, maxiter = 100, na.action = na.fail)
}

\arguments{
  \item{X}{a matrix or a data frame }
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{fixed.loc}{a logical, see details} 
  \item{location}{an optional vector giving the location of the data or the initial value for the location if it is estimated} 
  \item{init}{an optional starting value for the iteration}
  \item{steps}{fixed number of iteration steps to take, if \code{Inf} iteration is repeated until convergence (or until \code{maxiter} steps)}
  \item{eps}{tolerance for convergence}
  \item{maxiter}{maximum number of iteration steps. Ignored if \code{steps} 
 is finite} 
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
   } 

\details{\code{sign.shape} (Tyler's shape matrix), \code{symmsign.shape} (Duembgen's shape matrix), \code{rank.shape} and \code{signrank.shape} are the so called inner standardization matrices of location etc. tests based on spatial signs and ranks. When data is standardized using these matrices the corresponding sign or rank scores will appear \dQuote{uncorrelated}: the corresponding outer standardization matrices will be proportional to the identity matrix, see examples.

\code{spatial.shape} is a wrapper function for a unified access to all
four shape estimates. The choice of estimate is done via \code{score}:
\itemize{
\item \code{"sign"} for \code{signs.shape} 
\item \code{"symmsign"} for \code{symmsign.shape}
\item \code{"rank"} for \code{rank.shape} 
\item \code{"signrank"} for \code{signrank.shape}
}

\code{signrank.shape} and \code{sign.shape} include options to compute the shape matrix either with respect to fixed location (\code{fixed.loc = TRUE}) or so that the location and the shape are estimated simultaneously (\code{fixed.loc = FALSE}).

}

\value{
The estimate matrix with the (final estimate of or given) location vector 
as attribute \code{"location"}.
}

\author{Seija Sirkia, \email{seija.sirkia@iki.fi}, Jari Miettinen, \email{jari.p.miettinen@jyu.fi}}

\references{\cite{Oja, H., Randles R. (2004) Multivariate Nonparametric Tests. \emph{Statistical Science} 19, 598-605.}

\cite{Sirkia et al. (2009) Tests and estimates of shape based on spatial signs and ranks. \emph{Journal of Nonparametric Statistics}, 21, 155-176.}

\cite{Sirkia, S., Taskinen, S., Oja, H. (2007) Symmetrised M-estimators of scatter. \emph{Journal of Multivariate Analysis}, 98, 1611-1629.}}


\seealso{\code{\link[ICSNP]{tyler.shape}}, \code{\link[ICSNP]{duembgen.shape}}, also \link[=SCov]{spatial sign and rank covariance matrices} and \link[=spatial.rank]{spatial signs and ranks}}
\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rnorm(1500),ncol=3)\%*\%t(A)
symmsign.shape(X)
to.shape(symmsign.shape(X),trace=3)
spatial.shape(X,score="sign")
spatial.shape(X,score="sign",fixed.loc=TRUE)
to.shape(A\%*\%t(A))
# one-step shape estimate based on spatial ranks and covariance matrix:
spatial.shape(X,score="rank",init=cov(X),steps=1)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{htest}
