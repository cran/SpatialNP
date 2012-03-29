\name{Sphericity tests}
\alias{sr.sphere.test}
\title{Sphericity tests}
\description{
Tests of sphericity based on spatial signs and spatial signs of pairwise 
differences.
}

\usage{ sr.sphere.test(X, score = c("sign", "symmsign"), shape = NULL,
na.action = na.fail) }

\arguments{
  \item{X}{a matrix or a data frame}
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{shape}{a matrix with which the data should be standardized before the sphericity test}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}

\details{The test is for a null hypothesis of the form \dQuote{true
shape matrix is equal to the identity matrix}. Effectively, giving a
matrix as \code{shape} will produce a test of whether the true shape
is equal (in fact, proportional, since the scale of \code{shape} will
have no effect) to it. In that case the test will still be for
sphericity but the data is standardized beforehand.}

\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the statistic}
  \item{parameter}{the degrees of freedom for the statistic}
  \item{p.value}{the p-value for the test}
  \item{null.value}{the specified hypothesized value of the shape (always \code{diag(p)} where \code{p} is the number of dimensions)}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed}
  \item{data.name}{a character string giving the name of the data}
}


\author{Seija Sirkia, \email{seija.sirkia@iki.fi}}

\references{\cite{Sirkia et al. (2009) Tests and estimates of shape based on spatial signs and ranks. \emph{Journal of Nonparametric Statistics}, 21, 155-176.}}

\seealso{\link[=SCov]{sign and rank covariance matrices}}
\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rnorm(600),ncol=3)\%*\%t(A)
sr.sphere.test(X,score="sign")
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{htest}
