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

  t <- 1:(n - 1)

  # constructing CUSUM matrix
  rightmeans <- sweep(rightsums[, t, drop=FALSE], 2, n - t, '/')
  leftmeans <- sweep(leftsums[, t, drop=FALSE], 2, t, '/')
  cusum <- sweep(rightmeans - leftmeans, 2, sqrt(t * (n-t) / n), '*')
  return(cusum)
}

#' MissCUSUM transformation of a single vector with missing entries
#' @param x a vector with missing entries represented by NA
#' @return MissCUSUM transformed vector
#' @export
cusum.univariate.missing <- function(x){
  ob <- !is.na(x)
  z <- replace(x, !ob, 0)
  n <- length(x)
  leftsum <- cumsum(z)
  rightsum <- leftsum[n] - leftsum
  L <- cumsum(ob)
  R <- L[n] - L
  x.cusum <- (leftsum[-n] / L[-n] - rightsum[-n] / R[-n]) * sqrt(L[-n] * R[-n] / L[n])
  x.cusum[is.nan(x.cusum)] <- 0
  x.cusum
}

#' MissCUSUM transformation of a matrix with missing entries
#' @param x a matrix with missing entries represented by NA
#' @return MissCUSUM transformed matrix
#' @export
cusum.transform.missing <- function(x){
  if (!is.matrix(x)) {
    return(cusum.univariate.missing(x))
  } else if (ncol(x)==2){
    return(as.matrix(apply(x, 1, cusum.univariate.missing)))
  } else {
    return(t(apply(x, 1, cusum.univariate.missing)))
  }
}
