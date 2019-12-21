\name{rescale.variance}
\alias{rescale.variance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Noise standardisation for multivariate time series.
}
\description{
Each row of the input matrix is normalised by the estimated standard deviation computed through the median absolute deviation of increments.
}
\usage{
rescale.variance(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{
An input matrix of real values.
}
}
\details{
This is an auxiliary function used by the \code{InspectChangepoint} package.
}
\value{
A rescaled matrix of the same size is returned.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
x <- matrix(rnorm(40),5,8) * (1:5)
x.rescaled <- rescale.variance(x)
x.rescaled
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ts }% __ONLY ONE__ keyword per line
