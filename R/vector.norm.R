#' Norm of a vector
#' @description Calculate the entrywise L_q norm of a vector or a matrix
#' @param v a vector of real numbers
#' @param q a nonnegative real number or Inf
#' @param na.rm boolean, whether to remove NA before calculation
#' @return the entrywise L_q norm of a vector or a matrix
#' @export

vector.norm <- function(v, q=2, na.rm=FALSE){
    if (na.rm) v <- na.omit(v)
    if (q > 0) (sum(abs(v)^q))^(1/q)
    else if (q == 0) sum(v!=0)
    else if (q == Inf) max(abs(v))
    else NaN
}
