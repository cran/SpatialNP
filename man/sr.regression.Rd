\name{sr.regression}
\alias{sr.regression}
\title{Multivariate regression based on spatial signs or ranks}
\description{
This function fits linear models to data using spatial signs or ranks.
}

\usage{ sr.regression(formula, data=NULL, score = c("sign", "rank"),
ae = TRUE, eps = 1e-06, na.action = na.fail) }

\arguments{
  \item{formula}{a symbolic description of the model to be fit}
  \item{data}{an optional matrix or a data frame}
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{ae}{logical. Should the result be affine equivariant?}
  \item{eps}{tolerance for convergence}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}  } 

\details{ See \code{\link{formula}} or \code{\link{lm}} for more
details on specifying formulas.

The found matrix of coefficients \eqn{\beta}{\beta} for a model
\deqn{Y=X\beta^T+E}{Y=X \beta^T+E} is the one that solves \deqn{X^T
S(E(\beta)) = 0}{X^T S(E(\beta))} where \eqn{E(\beta)}{E(\beta)} is
the matrix of residuals corresponding to the choice of
\eqn{\beta}{\beta} and \eqn{S(E(\beta)}{S(E(\beta)} is the matrix of
the chosen scores (signs or ranks) of those residuals.

The coefficient estimates found based on ranks (\code{score = "rank"})
are always for the model with the intercept term included and the
intercept is computed separately. To emphasize this when rank based
coefficients are computed without the intercept term there is a row of
\code{NA} values in the result.}  

\value{ A matrix of estimated coefficients.  }


\author{Jaakko Nevalainen, \email{jaakko.nevalainen@uta.fi}}

\seealso{\code{\link[ICSNP]{tyler.shape}}, \code{\link{ae.hl.estimate}}}
\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rt(150,1),ncol=3)
Y<-X\%*\%t(A)+runif(150,-1,1)
sr.regression(Y~X,data.frame(Y=Y,X=X),score="sign")
sr.regression(Y~X,data.frame(Y=Y,X=X),score="rank")
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{regression}
