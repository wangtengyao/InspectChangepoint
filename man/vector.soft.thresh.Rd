\name{vector.soft.thresh}
\alias{vector.soft.thresh}
\title{
Soft thresholding operation
}
\description{
Performing soft thresholding to a vector or a matrix of real values, i.e. moving entries of the vector/matrix towards 0 by the threshold level.
}
\usage{
vector.soft.thresh(x, lambda)
}
\arguments{
  \item{x}{a vector or a matrix of real values}
  \item{lambda}{the threshold level, a non-negative real number}
}
\details{
This is an auxiliary function used by the \code{InspectChangepoint} package.
}
\value{
The thresholded vector/matrix is returned.
}
\examples{
n <- 10
lambda <- 1
vector.soft.thresh(rnorm(n), lambda)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~arith }% use one of  RShowDoc("KEYWORDS")
