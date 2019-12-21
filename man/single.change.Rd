\name{single.change}
\alias{single.change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High-dimensional time series with exactly one change in the mean structure
}
\description{
The data matrix is generated via X = mu + W, where mu is the mean structure matrix that captures the changepoint location and sparsity structure, and W is a random noise matrix.
}
\usage{
single.change(n, p, k, z, vartheta, sigma = 1, shape = 3, noise = 0, corr = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Time length of the observation
}
  \item{p}{
Dimension of the multivariate time series
}
  \item{k}{
Number of coordinates that undergo a change
}
  \item{z}{
Changepoint location, a number between 1 and n-1.
}
  \item{vartheta}{
The root mean squared change magnitude in coordinates that undergo a change
}
  \item{sigma}{
noise level
}
  \item{shape}{
How the signal strength is distributed across signal coordinates. When shape = 0, all signal coordinates are changed by the same amount; when shape = 1, their signal strength are proportional to 1, sqrt(2), ..., sqrt(k); when shape = 2, they are proportional to 1, 2, ..., k; when shape = 3, they are proportional to 1, 1/sqrt(2), ..., 1/sqrt(k).
}
  \item{noise, corr}{
Noise structure of the multivarite time series. For noise = 0, 0.5, 1, columns of W have independent multivariate normal distribution with covariance matrix Sigma. When noise = 0, Sigma = sigma^2 * I_p; when noise = 0.5, noise has local dependence structure given by Sigma_{i,j} = sigma*corr^|i-j|; when noise = 1, noise has global dependence structure given by matrix(corr,p,p)+diag(p)*(1-corr))) * sigma. When noise = 2, rows of the W are independent and each having an AR(1) structure given by W_{j,t} = W_{j,t-1} * sqrt(corr) + rnorm(sd = sigma) * sqrt(1-corr). For noise = 3, 4, entries of W have i.i.d. uniform distribution and exponential distribution respectively, each centred and rescaled to have zero mean and variance sigma^2.
}
}
\value{
An S3 object of the class 'hdchangeseq' is returned.
  \item{x }{The generated data matrix}
  \item{mu }{The mean structure of the data matrix}
}

\seealso{
\code{\link{plot.hdchangeseq}}
}
\examples{
n <- 2000; p <- 100; k <- 10; z <- 800; vartheta <- 1; sigma <- 1; shape <- 3
noise <- 0; corr <- 0
obj <- single.change(n,p,k,z,vartheta,sigma, shape, noise, corr)
plot(obj, noise = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~data }% use one of  RShowDoc("KEYWORDS")
