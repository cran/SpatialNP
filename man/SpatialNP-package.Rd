\name{SpatialNP-package}
\alias{SpatialNP-package}
\alias{SpatialNP}
\docType{package}
\title{
Multivariate nonparametric tests and estimates based on spatial signs and ranks.
}
\description{
This package contains tests and an estimate of location, tests of independence,
tests of sphericity and several estimates of shape all based on spatial signs, 
symmetrized signs (that is, signs of pairwise differences), ranks and signed ranks.
}
\details{
\tabular{ll}{
Package: \tab SpatialNP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-07-02\cr
License: \tab GPL (>= 2)\cr
}

In this package there are three functions for inference,
\code{\link{sr.loc.test}}, \code{\link{sr.indep.test}} and
\code{\link{sr.sphere.test}}, for location, independence and
sphericity tests. The so called \link[=rank.shape]{inner} and
\link[=SCov]{outer} standardization matrices are also available as
well as the actual \link[=spatial.rank]{sign and rank} score
functions, together with a utility function \code{\link{to.shape}}.  }


\author{
Seija Sirkia, Klaus Nordhausen, Hannu Oja
Maintainer: Seija Sirkia, \email{seija.sirkia@iki.fi}}


\references{\cite{Oja, H., Randles R. (2004) Multivariate Nonparametric Tests. \emph{Statistical Science} 19, 598-605}}
\keyword{ package }
\seealso{
package \pkg{\link[ICSNP:ICSNP-package]{ICSNP}}}


