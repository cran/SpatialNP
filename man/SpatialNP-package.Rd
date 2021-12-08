\name{SpatialNP-package}
\alias{SpatialNP-package}
\alias{SpatialNP}
\docType{package}
\title{
Multivariate nonparametric tests and estimates based on spatial signs and ranks.
}
\description{
Test and estimates of location, tests of independence, tests of sphericity and several estimates of shape all based on spatial signs, symmetrized signs, ranks and signed ranks. For details, see Oja and Randles (2004) <doi:10.1214/088342304000000558> and Oja (2010)  <doi:10.1007/978-1-4419-0468-3>.
}
\details{
\tabular{ll}{
Package: \tab SpatialNP\cr
Type: \tab Package\cr
Version: \tab 1.1-5\cr
Date: \tab 2021-12-08\cr
License: \tab GPL (>= 2)\cr
}

There are three functions for inference,
\code{\link{sr.loc.test}}, \code{\link{sr.indep.test}} and
\code{\link{sr.sphere.test}}, for location, independence and
sphericity tests. The so called \link[=rank.shape]{inner} and
\link[=SCov]{outer} standardization matrices are also available as
well as the actual \link[=spatial.rank]{sign and rank} score
functions, together with a utility function \code{\link{to.shape}}.  }


\author{
Seija Sirkia, Jari Miettinen, Klaus Nordhausen, Hannu Oja, Sara Taskinen \cr
Maintainer: Klaus Nordhausen, \email{klaus.k.nordhausen@jyu.fi}}

\keyword{ package }
\seealso{
package \pkg{\link[ICSNP:ICSNP-package]{ICSNP}} \cr
package \pkg{\link[MNM:MNM-package]{MNM}}}


