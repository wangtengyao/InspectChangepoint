#' Power method for finding the leading left singular vector of a matrix
#' @param A a rectangular
#' @param eps tolerance for convergence (in Frobenius norm)
#' @param maxiter maximum iteration
#' @return a unit-length leading eigenvector of A
#' @export
power.method <- function(A, eps = 1e-10, maxiter = 10000){
  if (nrow(A) != ncol(A)) stop('powerMethod requires a square matrix')
  d <- nrow(A)
  v <- rnorm(d); v <- v/vector.norm(v)

  for (i in seq_len(maxiter)){
    v_new <- A%*%v; v_new <- v_new/vector.norm(v_new)
    if ((vector.norm(v_new - v)) < eps) break
    v <- v_new
  }
  if (i == maxiter) warning('max iter reached without convergence.')
  return(v_new)
}

