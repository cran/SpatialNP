\name{mvhuberM}
\alias{mvhuberM}
\alias{symmhuber}

\title{Multivariate Huber's M-estimator and its symmetrized version}

\description{
Iterative algorithms to estimate M-estimators of location and scatter as well as symmetrized M-estimator using Huber's weight functions.

}
\usage{
mvhuberM(X, qg = 0.9, fixed.loc = FALSE, location = NULL, init =
NULL, steps = Inf, eps = 1e-06, maxiter = 100, na.action = na.fail)

symmhuber(X, qg = 0.9, init = NULL, steps = Inf, eps = 1e-6, 
maxiter = 100, na.action = na.fail)
}

\arguments{
  \item{X}{a matrix or a data frame}
  \item{qg}{a tuning parameter. The default is 0.9, see details}
  \item{fixed.loc}{a logical, see details}
  \item{location}{an optional vector giving the location of the data or the initial value for the location if it is estimated}
  \item{init}{an optional starting value for scatter}  
  \item{steps}{fixed number of iteration steps to take, if \code{Inf} iteration is repeated until convergence (or until \code{maxiter} steps)}
  \item{eps}{tolerance for convergence}
  \item{maxiter}{maximum number of iteration steps. Ignored if \code{steps} 
 is finite} 
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail.}
}
\details{
\code{mvhuberM} computes multivariate M-estimators of location and scatter 
using Huber's weight functions. The tuning parameter \code{qg} defines cutoff-point c for weight functions so that \eqn{c=F^{-1}(q)}{c=F^{-1}(q)}, where F is the cdf of \eqn{\chi^2}{chi^2}-distribution with p degrees of freedom. The estimators with maximal breakdown point are obtained with the choice qg=F(p+1). If \code{fixed.loc} is set TRUE, scatter estimator is computed with fixed location given by 
\code{location} (default is column means). 

\code{symmhuber} computes Huber's M-estimator of scatter using pairwise
differences of the data therefore avoiding location estimation.   
}

\value{
\code{mvhuberM} returns a list with components
\item{location }{a vector}
\item{scatter }{a matrix}

\code{symmhuber} returns a matrix.
}

\references{
\cite{Huber, P.J. (1981), Robust Statistics, Wiley, New York.}

\cite{Lopuhaa, H.P. (1989). On the relation between S-estimators and M-estimators of multivariate location and covariance. \emph{Annals of Statistics}, 17, 1662-1683.}

\cite{Sirkia, S., Taskinen, S., Oja, H. (2007) Symmetrised M-estimators of scatter. \emph{Journal of Multivariate Analysis}, 98, 1611-1629.}
}
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}, \cr
Jari Miettinen, \email{jari.p.miettinen@jyu.fi}} 

\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rnorm(1500),ncol=3)\%*\%t(A)
mvhuberM(X)
symmhuber(X)
}
\keyword{multivariate}
\keyword{robust}
