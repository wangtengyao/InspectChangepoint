\name{multi.change}
\alias{multi.change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High-dimensional time series with multiple changepoints
}
\description{
The data matrix is generated via X = mu + W, where mu is the mean structure matrix that captures the changepoint locations and sparsity structure, and W is a random noise matrix having independent N(0,sigma^2) entries. 
}
\usage{
multi.change(n, p, ks, zs, varthetas, sigma = 1, overlap = 0, shape = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Time length of the observation
}
  \item{p}{
Dimension of the multivariate time series
}
  \item{ks}{
A vector describing the number of coordinates that undergo a change in each changepoint. If only a scalar is supplied, each changepoint will have the same number of coordinates that undergo a change.
}
  \item{zs}{
A vector describing the locations of the changepoints.
}
  \item{varthetas}{
A vector describing the root mean squared change magnitude in coordinates that undergo a change for each changepoint. If only a scalar is supplied, each changepoint will have the same signal strength value. 
}
  \item{sigma}{
noise level
}
  \item{overlap}{
A number between 0 and 1. The proportion of overlap in the signal coordinates for successive changepoints. 
}
  \item{shape}{
How the signal strength is distributed across signal coordinates. When shape = 0, all signal coordinates are changed by the same amount; when shape = 1, their signal strength are proportional to 1, sqrt(2), ..., sqrt(k); when shape = 2, they are proportional to 1, 2, ..., k; when shape = 3, they are proportional to 1, 1/sqrt(2), ..., 1/sqrt(k).
}
}
\value{
An S3 object of the class 'hdchangeseq' is returned.
  \item{x }{The generated data matrix}
  \item{mu }{The mean structure of the data matrix}
%% ...
}
\seealso{
\code{\link{plot.hdchangeseq}}
}
\examples{
n <- 2000; p <- 200; ks <- 40; zs <- c(500,1000,1500); varthetas <- c(0.1,0.15,0.2); overlap <- 0.5
obj <- multi.change(n, p, ks, zs, varthetas, overlap)
plot(obj, noise = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~data }% use one of  RShowDoc("KEYWORDS")
