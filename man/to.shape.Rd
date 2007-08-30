\name{to.shape}
\alias{to.shape}
\title{Rescale a matrix to a shape matrix}
\description{
This function rescales a given matrix such that its determinant, trace or the value of the first diagonal element meets a criteria.}
}
\usage{
to.shape(M, determ, trace, first)
}
\arguments{
  \item{M}{a matrix to be scaled}
  \item{determ}{required value for determinant}
  \item{trace}{required value for trace}
  \item{first}{required value of the first diagonal element}
}
}

\details{ If \code{determ}, \code{trace} or \code{first} is given
\code{M} is scaled such that its determinant, trace or first diagonal
element, respectively, equals that value. If none of the three is
given \code{M} is scaled such that its determinant equals one. If more
than one criteria is given the first of them is used and the others
silently ignored.  }

\value{
The rescaled matrix
}

\references{Paindaveine D. (2007) A Canonical Definition of Shape. Statistics and Probability Letters (tentatively accpeted)}
\author{Seija Sirkia, \email{ssirkia@maths.jyu.fi}}

%\author{ ~~who you are~~ }
\note{A shape matrix is symmetric and positive definite square matrix. In order for the result to be such the argument matrix \code{M} should also be symmetric and positive definite square matrix. However, the function does not check for this. Expect to see errors if \code{M} is of inappropriate type.}

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
(A<-matrix(c(1,2,-3,4,3,-2,-1,0,4),ncol=3))
to.shape(A\%*\%t(A))
to.shape(A\%*\%t(A),trace=3)
to.shape(A\%*\%t(A),first=1)

}
\keyword{multivariate}
