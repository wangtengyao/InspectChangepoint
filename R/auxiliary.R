###########   Matrices and vectors   ##############

#' Norm of a vector
#' @description Calculate the entrywise L_q norm of a vector or a matrix
#' @param v a vector of real numbers
#' @param q a nonnegative real number or Inf
#' @param na.rm boolean, whether to remove NA before calculation
#' @return the entrywise L_q norm of a vector or a matrix
vector.norm <- function(v, q = 2, na.rm = FALSE){
  if (na.rm) v <- na.omit(v)
  M <- max(abs(v))
  if (M == 0) return(0) else v <- v/M
  if (q == Inf) {
    nm <- max(abs(v))
  } else if (q > 0) {
    nm <- (sum(abs(v)^q))^(1/q)
  } else if (q == 0) {
    nm <- sum(v!=0)
  } else {
    return(NaN)
  }
  return(nm * M)
}

#' Normalise a vector
#' @param v a vector of real numbers
#' @param q a nonnegative real number or Inf
#' @param na.rm boolean, whether to remove NA before calculation
#' @return normalised version of this vector
vector.normalise <- function(v, q = 2, na.rm = FALSE){
  return(v / vector.norm(v, q, na.rm))
}

#' Clipping a vector from above and below
#' @description Clipping vector or matrix x from above and below
#' @param x a vector of real numbers
#' @param upper clip above this value
#' @param lower clip below this value
#' @return the entrywise L_q norm of a vector or a matrix
vector.clip <- function(x, upper = Inf, lower = -upper){
  if (upper < lower)  stop("upper limit cannot be below lower limit")
  x[x<lower]<-lower;
  x[x>upper]<-upper;
  x
}

#' Soft thresholding a vector
#' @param x a vector of real numbers
#' @param lambda soft thresholding value
#' @return a vector of the same length
#' @description entries of v are moved towards 0 by the amount lambda until they hit 0.
vector.soft.thresh <- function(x, lambda){
  sign(x)*pmax(0,(abs(x)-lambda))
}

#' Generate a random unit vectors in R^n
#' @param n length of random vector
random.UnitVector <- function(n){
  v = rnorm(n)
  v/vector.norm(v)
}

#' Noise standardisation for multivariate time series.
#' @description Each row of the input matrix is normalised by the estimated standard deviation computed through the median absolute deviation of increments.
#' @param x An input matrix of real values.
#' @details This is an auxiliary function used by the \code{InspectChangepoint} package.
#' @return A rescaled matrix of the same size is returned.
#' @examples
#' x <- matrix(rnorm(40),5,8) * (1:5)
#' x.rescaled <- rescale.variance(x)
#' x.rescaled

rescale.variance <- function(x){
  p <- dim(x)[1]
  n <- dim(x)[2]
  for (j in 1:p){
    v <- x[j,]
    v <- v[!is.na(v)]
    scale <- mad(diff(v))/sqrt(2)
    x[j,] <- x[j,] / scale
  }
  return(x)
}

#' Print percentage
#' @param ind a vector of for loop interator
#' @param tot a vector of for loop lengths
#' @return on screen output of percentage
printPercentage <- function (ind, tot){
  ind <- as.vector(ind); tot <- as.vector(tot)
  if ((length(tot) > 1) & (length(ind) == 1)) {ind <- match(ind, tot); tot <- length(tot)}
  len <- length(ind)
  contrib <- rep(1,len)
  if (len > 1) {
    for (i in (len-1):1) contrib[i] <- contrib[i+1] * tot[i+1]
  }
  grand_tot <- contrib[1] * tot[1]
  count <- (sum(contrib * (ind - 1)) + 1)
  out <- ""
  if (sum(ind-1)>0) out <- paste0(rep("\b", nchar(round((count-1)/grand_tot * 100))+1), collapse = "")
  out <- paste0(out, round(count/grand_tot*100), "%")
  if (identical(ind, tot)) out <- paste0(out, '\n')
  cat(out)
  return(NULL)
}
