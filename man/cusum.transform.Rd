\name{cusum.transform}
\alias{cusum.transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CUSUM transformation
}
\description{
Performing CUSUM transformation to the input matrix of multivariate time series. If the input is a vector, it is treated as a matrix of one row. 
}
\usage{
cusum.transform(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{
Input matrix or vector.
}
}
\details{
For any integers p and n, the CUSUM transformation \eqn{T_{p,n}: R^{p\times n}\to R^{p\times (n-1)}} is defined by
\deqn{
[T_{p,n}(M)]_{j,t} := \sqrt{t(n-t)/n}\biggl(\frac{1}{n-t}\sum_{r=t+1}^n M_{j,r} - \frac{1}{t}\sum_{r=1}^t M_{j,r}\biggr).
}
}
\value{
The transformed matrix is returned. Note that the returned matrix has the same number of rows but one fewer columns compared with the input matrix.
}
\examples{
x <- matrix(rnorm(20),4,5)
cusum.transform(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ts }% use one of  RShowDoc("KEYWORDS")
