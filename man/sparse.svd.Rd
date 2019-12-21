\name{sparse.svd}
\alias{sparse.svd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse left leading singular vector
}
\description{
Estimating the sparse left leading singular vector by first computing a maximiser Mhat of the convex problem
\deqn{<Z, M> - \lambda |M|_1} 
subject to the Schatten norm constraint |M|_schatten <= 1 using alternating direction method of multipliers (ADMM). Then the leading left singular vector of Mhat is returned.
}
\usage{
sparse.svd(Z, lambda, schatten = c(1, 2), tolerance = 1e-05, max.iter = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z }{
Input matrix whose left leading singular vector is to be estimated.
}
  \item{lambda }{
Regularisation parameter
}
  \item{schatten }{
Schatten norm constraint to be used. Default uses Schatten-2-norm, i.e. the Frobenius norm. Also possible to use Schatten-1-norm, the nuclear norm.
}
  \item{tolerance }{
Tolerance criterion for convergence of the ADMM algorithm.
}
  \item{max.iter }{
Maximum number of iteration in the ADMM algorithm.
}
}
\details{
In case of schatten = 2, a closed-form solution for Mhat using matrix soft thresholding is possible. We use the closed-form solution instead of the ADMM algorithm to speed up the computation.
}
\value{
A vector that has the same length as number of rows of the input matrix is returned.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
Z <- matrix(rnorm(20),4,5)
lambda <- 0.5
sparse.svd(Z, lambda)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~algebra }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~optimize }% __ONLY ONE__ keyword per line
