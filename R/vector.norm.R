#' Norm of a vector
#' @description Calculate the entrywise L_q norm of a vector or a matrix
#' @param v a vector of real numbers
#' @param q a nonnegative real number or Inf
#' @param na.rm boolean, whether to remove NA before calculation
#' @return the entrywise L_q norm of a vector or a matrix
#' @export

vector.norm <- function (v, q = 2, na.rm = FALSE)
{
  if (na.rm)
    v <- na.omit(v)
  M <- max(abs(v))
  if (M == 0)
    return(0)
  else v <- v/M
  if (q == Inf) {
    nm <- max(abs(v))
  }
  else if (q > 0) {
    nm <- (sum(abs(v)^q))^(1/q)
  }
  else if (q == 0) {
    nm <- sum(v != 0)
  }
  else {
    return(NaN)
  }
  return(nm * M)
}
