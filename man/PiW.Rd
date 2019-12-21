\name{PiW}
\alias{PiW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Projection onto the standard simplex
}
\description{
The input vector is projected onto the standard simplex, i.e. the set of vectors of the same length as the input vector with non-negative entries that sum to 1.
}
\usage{
PiW(v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v }{
Input vector
}
}
\details{
This is an auxiliary function used by the \code{InspectChangepoint} package.
}
\value{
A vector in the standard simplex that is closest to the input vector is returned.
}
\references{
Chen, Y. and Ye, X. (2011) Projection onto a simplex. arXiv preprint, arxiv:1101.6081.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
v <- rnorm(10)
PiW(v)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~algebra }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~optimize }% __ONLY ONE__ keyword per line
