#' CUSUM transformation
#' @description Performing CUSUM transformation to the input matrix of multivariate time series. If the input is a vector, it is treated as a matrix of one row.
#' @details For any integers p and n, the CUSUM transformation \eqn{T_{p,n}: R^{p\times n}\to R^{p\times (n-1)}} is defined by
#' \deqn{
#'    [T_{p,n}(M)]_{j,t} := \sqrt{t(n-t)/n}\biggl(\frac{1}{n-t}\sum_{r=t+1}^n M_{j,r} - \frac{1}{t}\sum_{r=1}^t M_{j,r}\biggr).
#' }
#'
#' @param x input matrix
#' @return The transformed matrix is returned. Note that the returned matrix has the same number of rows but one fewer columns compared with the input matrix.
#' @examples
#' x <- matrix(rnorm(20),4,5)
#' cusum.transform(x)
#' @export

cusum.transform <- function(x){
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation

    leftsums <- t(apply(x, 1, cumsum))
    rightsums <- leftsums[, n] - leftsums
    leftsums <- t(leftsums)
    rightsums <- t(rightsums)
    t <- 1:(n - 1)

    # constructing CUSUM matrix
    return(t((rightsums[t,] / (n-t) - leftsums[t,] / t) * sqrt(t * (n-t) / n)))
}
