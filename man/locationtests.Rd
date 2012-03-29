\name{Location tests}
\alias{sr.loc.test}
\title{Spatial sign and rank tests of multivariate location}

\description{ Multivariate tests of location of one or more samples
based on spatial signs and (signed) ranks. In case of one sample the
null hypothesis about a given location is tested. In case of several
samples the null hypothesis is that all samples have the same
location.  }

\usage{ sr.loc.test(X, Y = NULL, g = NULL, score = c("sign", "rank"),
nullvalue = NULL, cond = FALSE, cond.n = 1000,
na.action = na.fail,...)  }

\arguments{
  \item{X}{a matrix or a data frame }
  \item{Y}{an optional matrix or a data frame }
  \item{g}{a factor giving the groups (may contain just one level) }
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{nullvalue}{location to be tested in the one sample case (ignored if there is more than one sample)}
  \item{cond}{logical. Should the conditionally distribution free test be used? (Ignored if \code{score} is \code{"rank"})}
  \item{cond.n}{number of permutations to use in the conditionally distribution free test}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
  \item{\dots}{further arguments to be passed to other functions}
}
\details{
\code{X} should contain the the whole data set and \code{g} should describe the groups, or, if there is only one group, \code{g} may be missing. Alternatively, if there are two samples \code{X} may contain only the first sample while the second sample is given in \code{Y} and \code{g} is ignored. Note that in the one sample case when \code{rank} is chosen as \code{score} the function in fact uses signed ranks.

Note that the conditionally distribution free p-value is only provided for the sign based version of the test.
}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the statistic}
  \item{parameter}{the degrees of freedom for the statistic or the number of replications if conditionally distribution free p-value was used}
  \item{p.value}{the p-value for the test}
  \item{null.value}{the specified hypothesized value of the (common) location}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed}
  \item{data.name}{a character string giving the name of the data (and grouping vector)}
}

\author{Seija Sirkia, \email{seija.sirkia@iki.fi}}

\references{\cite{Oja, H., Randles R. (2004) Multivariate Nonparametric Tests. \emph{Statistical Science} 19, 598-605.}}

\seealso{\link[=spatial.rank]{Spatial signs and ranks}}

\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-rbind(matrix(rnorm(1500),ncol=3),matrix(rnorm(750)+1,ncol=3))\%*\%t(A)
sr.loc.test(X,cond=TRUE)
X[1:250,]<-X[1:250,]+1
g<-factor(rep(c(1,2,3),each=250))
sr.loc.test(X,g=g,score="rank")
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{htest}
