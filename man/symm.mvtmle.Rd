\name{symm.mvtmle}
\alias{symm.mvtmle}
\alias{symm.mvtmle.inc}

\title{Symmetrized M-estimators of scatter with the weights of the t-distribution}

\description{
Iterative algorithms to estimate symmetrized M-estimators of scatter using weights of the t-distribution.

}
\usage{
symm.mvtmle(X, nu=1, init=NULL, steps=Inf, eps=1e-6, 
maxiter=100, na.action = na.fail)

symm.mvtmle.inc(X, nu=1, m=10, init=NULL, steps=Inf, permute=TRUE,
eps=1e-6, maxiter=100, na.action = na.fail)
}

\arguments{
  \item{X}{a matrix or a data frame}
  \item{nu}{the degrees of freedom of the t-distribution. The default is 1. Must be larger than 0.}
  \item{init}{an optional starting value for scatter}  
  \item{steps}{fixed number of iteration steps to take, if \code{Inf} iteration is repeated until convergence (or until \code{maxiter} steps)}
  \item{m}{a parameter in \code{symm.mvtmle.inc} which defines how many pairwise differences are used, see details.}
  \item{permute}{logical in \code{symm.mvtmle.inc} which determines whether the rows of \code{X} are permuted randomly, see details.}
  \item{eps}{tolerance for convergence}
  \item{maxiter}{maximum number of iteration steps. Ignored if \code{steps} 
 is finite} 
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail.}
}
\details{
\code{symm.mvtmle} computes M-estimator of scatter using weights of the t-distribution and pairwise
differences of the data. Hence, location estimation is not needed.   

\code{symm.mvtmle.inc} is a computationally lighter estimator to approximate symmetrized M-estimator of scatter which uses weights of the t-distribution. 
Only a subset of the pairwise differences are used in the computation in the incomplete case. The magnitude of the subset used is controlled
by the argument \code{m} which is half of the number of how many differences each observation is part of. Differences of successive observations are used, 
and therefore random permutation of the rows of \code{X} is suggested and is the default choice in the function. For details see Miettinen et al., 2016.
}

\value{
\code{symm.mvtmle} returns a matrix.

\code{symm.mvtmle.inc} returns a matrix.
}

\references{
\cite{Huber, P.J. (1981), Robust Statistics, Wiley, New York.}

\cite{Sirkia, S., Taskinen, S., Oja, H. (2007) Symmetrised M-estimators of scatter. \emph{Journal of Multivariate Analysis}, 98, 1611-1629.}

\cite{Duembgen, L., Pauly, M., Schweizer, T. (2015) M-Functionals of multivariate scatter.
\emph{Statistics Surveys 9}, 32-105.}

\cite{Miettinen, J., Nordhausen, K., Taskinen, S., Tyler, D.E. (2016) On the computation of symmetrized M-estimators of scatter. \emph{In Agostinelli, C. Basu, A., Filzmoser, P. and Mukherje, D. (editors) ''Recent Advances in Robust Statistics: Theory and Application''}, 131-149, Springer India, New Delhi.}
}
\author{Jari Miettinen, \email{jari.p.miettinen@aalto.fi}, \cr
Klaus Nordhausen, \email{klaus.nordhausen@tuwien.ac.at}
} 

\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rnorm(1500),ncol=3)\%*\%t(A)
symm.mvtmle(X, nu=2)
symm.mvtmle.inc(X, nu=2, m=20)
}
\keyword{multivariate}
\keyword{robust}
