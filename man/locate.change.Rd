\name{locate.change}
\alias{locate.change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single changepoint estimation
}
\description{
Estimate the location of one changepoint in a multivariate time series. It uses the function \code{\link{sparse.svd}} to estimate the best projection direction, then using univariate CUSUM statistics of the projected time series to estimate the changepoint location.
}
\usage{
locate.change(x, lambda, schatten = 2, sample.splitting = FALSE, 
standardize.series = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix of multivariate time series
}
  \item{lambda}{
Regularisation parameter. If no value is supplied, the dafault value is chosen to be log(log(n)*p/2) for p and n number of rows and columns of the data matrix x respectively.
}
  \item{schatten}{
The Schatten norm constraint to use in the \code{\link{sparse.svd}} function. Default is schatten = 2, i.e. a Frobenius norm constraint.
}
  \item{sample.splitting}{
  Whether the changepoint should be estimated via sample splitting. The theoretical result is proven only for the sample splitted version of the algorithm. However, the default setting in practice is without sample splitting.
  }
  \item{standardize.series}{
Whether the given time series should be standardised before estimating the projection direction. Default is FALSE, i.e. the input series is assume to have variance 1 in each coordinate.
}
}
\value{
\item{changepoint}{A single integer value estimate of the changepoint location is returned. If the estimated changepoint is z, it means that the multivariate time series is piecewise constant up to z and from z+1 onwards.}
\item{cusum}{The maximum absolute CUSUM statistic of the projected univariate time series associated with the estimated changepoint.}
}
\references{
Wang, T., Samworth, R. J. (2016) High-dimensional changepoint estimation via sparse projection. Arxiv preprint: arxiv1606.06246.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
n <- 2000; p <- 1000; k <- 32; z <- 400; vartheta <- 0.12; sigma <- 1; shape <- 3
noise <- 0; corr <- 0
obj <- single.change(n,p,k,z,vartheta,sigma,shape,noise,corr)
x <- obj$x
locate.change(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ts }% use one of  RShowDoc("KEYWORDS")

