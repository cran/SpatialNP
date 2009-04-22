\name{Independence tests}
\alias{sr.indep.test}

\title{Multivariate test of independence based on spatial signs or
ranks} 

\description{ Test of independence between two sets of
variables. Inference is based on the spatial signs of the
observations, symmetrized signs of the observations or spatial signed
ranks of the observations.  }

\usage{ sr.indep.test(X, Y = NULL, g = NULL, score = c("sign",
"symmsign", "rank"), regexp= FALSE, cond = FALSE, cond.n = 1000,
na.action = na.fail) }

\arguments{
  \item{X}{a matrix or a data frame}
  \item{Y}{an optional matrix or a data frame}
  \item{g}{a factor giving the two sets of variables, or numeric vector or vector of column names giving the first set of variables. See details}
  \item{score}{a character string indicating which transformation of the observations should be used}
  \item{regexp}{logical. Is \code{g} a regular expression?}
  \item{cond}{logical. Should the conditionally distribution free test be used?}
  \item{cond.n}{Number of permutations to use in the conditionally distribution free test}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}

\details{\code{X} should contain the first set of variables and
\code{Y} the second with matching rows. Alternatively, \code{X} should
contain both sets and \code{g} should be a factor of length equal to
number of columns of \code{X}, or, \code{g} should be a numeric or
character vector naming the variables in the first set. If \code{g} is
a character vector it is assumed to name all wanted columns exactly,
unless \code{regexp} is \code{TRUE}.}

\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the statistic}
  \item{parameter}{the degrees of freedom for the statistic or the number of replications if conditionally distribution free p-value was used}
  \item{p.value}{the p-value for the test}
  \item{null.value}{the specified hypothesized value of the measure of dependence (always 0)}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed}
  \item{data.name}{a character string giving the name of the data (and grouping vector)}
}

\author{Seija Sirkia, \email{seija.sirkia@iki.fi}}

\seealso{\link[=spatial.rank]{Spatial signs and ranks}}
\examples{
A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3)
X<-matrix(rt(150,1),ncol=3)\%*\%t(A)
Y<-cbind(X+runif(150,-1,1),runif(50))
sr.indep.test(X,Y)
#alternative calls:
Z<-cbind(X,Y)
colnames(Z)<-c("a1","a2","a3","b1","b2","b3","b4")
g<-factor(c(rep(1,3),rep(2,4)))
sr.indep.test(Z,g=g)
sr.indep.test(Z,g=c("b"),regexp=TRUE)
sr.indep.test(Z,g=1:3)

}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{htest}
