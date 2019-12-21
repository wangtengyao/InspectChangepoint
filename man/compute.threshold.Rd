\name{compute.threshold}
\alias{compute.threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing threshold used in \code{inspect}
}
\description{
The threshold level to be used in \code{inspect} is computed via Monte Carlo simulation of multivariate time series that do not contain any changepoints.
}
\usage{
compute.threshold(n, p, nrep = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Time length of the observation.
}
  \item{p}{
Dimension of the multivariate time series.
}
  \item{nrep}{
Number of Monte Carlo repetition to be used.
}
}
\value{
A numeric value indicating the threshold level that should be used based on the Monte Carlo simulation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
compute.threshold(n = 200, p = 50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ts }% use one of  RShowDoc("KEYWORDS")
