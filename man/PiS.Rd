\name{PiS}
\alias{PiS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix projection onto the nuclear norm unit sphere
}
\description{
Projection (with respect to the inner product defined by the Frobenius norm) of a matrix onto the unit sphere defined by the nuclear norm.
}
\usage{
PiS(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M }{
Input matrix
}
}
\details{
This is an auxiliary function used by the \code{InspectChangepoint} package. The projection is achieved by first performing a singular value decomposition, then projecting the vector of singular values onto the standard simplex, and finally using singular value decomposition in reverse to build the projected matrix.
}
\value{
A matrix of the same dimension as the input is returned.
}

\examples{
M <- matrix(rnorm(20),4,5)
PiS(M)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~algebra }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~optimize }% __ONLY ONE__ keyword per line
