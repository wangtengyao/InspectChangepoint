\name{vector.norm}
\alias{vector.norm}
\title{
Vector norms and entrywise matrix norms
}
\description{
Compute the (entrywise) L_q norm of a vector or a matrix
}
\usage{
vector.norm(v, q = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
input vector or matrix
}
  \item{q}{
L_q norm to use
}
}
\details{
This is an auxiliary function used by the \code{InspectChangepoint} package.
}
\value{
A non-negative numeric value is returned for the L_q norm of v. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
v <- rnorm(10)
vector.norm(v,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~arith }% use one of  RShowDoc("KEYWORDS")
