#' Single changepoint estimation
#' @description Estimate the location of one changepoint in a multivariate time
#' series. It uses the function \code{\link{sparse.svd}} to estimate the best
#' projection direction, then using univariate CUSUM statistics of the projected
#' time series to estimate the changepoint location.
#' @param x A (p x n) data matrix of multivariate time series, each column
#' represents a data point
#' @param lambda Regularisation parameter. If no value is supplied, the dafault
#' value is chosen to be sqrt(log(log(n)*p/2)) for p and n number of rows and
#' columns of the data matrix x respectively.
#' @param schatten The Schatten norm constraint to use in the \code{\link{sparse.svd}}
#'  function. Default is schatten = 2, i.e. a Frobenius norm constraint.
#' @param sample.splitting Whether the changepoint should be estimated via
#' sample splitting. The theoretical result is proven only for the sample
#' splitted version of the algorithm. However, the default setting in practice
#' is without sample splitting.
#' @param standardize.series Whether the given time series should be
#' standardised before estimating the projection direction. Default is FALSE,
#' i.e. the input series is assume to have variance 1 in each coordinate.
#' @param view.cusum Whether to show a plot of the projected CUSUM series
#' @return A list of two items:
#' \itemize{
#'   \item changepoint - A single integer value estimate of the changepoint
#'   location is returned. If the estimated changepoint is z, it means that the
#'   multivariate time series is piecewise constant up to z and from z+1
#'   onwards.
#'   \item cusum - The maximum absolute CUSUM statistic of the projected
#'   univariate time series associated with the estimated changepoint.
#'   \item vector.proj - the vector of projection, which is proportional to an estimate of the vector of change.
#' }
#' @references Wang, T., Samworth, R. J. (2016) High-dimensional changepoint estimation via sparse projection. Arxiv preprint: arxiv1606.06246.
#' @examples
#' n <- 2000; p <- 1000; k <- 32; z <- 400; vartheta <- 0.12; sigma <- 1; shape <- 3
#' noise <- 0; corr <- 0
#' obj <- single.change(n,p,k,z,vartheta,sigma,shape,noise,corr)
#' x <- obj$x
#' locate.change(x)
#' @export

locate.change <- function(x, lambda, schatten=2, sample.splitting=FALSE,
                          standardize.series=FALSE, view.cusum=FALSE)
{
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation
    if (missing(lambda)) lambda <- sqrt(log(log(n)*p)/2)
    if (standardize.series) x <- rescale.variance(x)
    if (sample.splitting){
        x1 <- x[,seq(1,n,by=2)]
        x2 <- x[,seq(2,n,by=2)]
    } else {
        x1 <- x
        x2 <- x
    }

    # construct cusum matrix of x
    cusum.matrix1 <- cusum.transform(x1)
    if (sample.splitting) {
        cusum.matrix2 <- cusum.transform(x2)
    } else {
        cusum.matrix2 <- cusum.matrix1
    }

    # estimate changepoint
    if (lambda >= max(abs(cusum.matrix1))) lambda <- max(abs(cusum.matrix1)) - 1e-10

    vector.proj <- sparse.svd(cusum.matrix1, lambda, schatten);
    cusum.proj <- t(cusum.matrix2)%*%vector.proj

    if (view.cusum) plot(as.numeric(cusum.proj), ylab='projected cusum', pch=20)

    ret <- NULL
    ret$changepoint <- which.max(abs(cusum.proj))
    if (sample.splitting) ret$changepoint <- ret$changepoint * 2
    ret$cusum <- max(abs(cusum.proj))
    ret$vector.proj <- vector.proj

    return(ret)
}
