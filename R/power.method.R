#' Power method for finding the leading eigenvector of a symmetric matrix
#' @param A a square symmetric matrix
#' @param eps tolerance for convergence (in Frobenius norm)
#' @param maxiter maximum iteration
#' @return a unit-length leading eigenvector of A
#' @export

power.method <- function(A, eps = 1e-10, maxiter = 10000){
  if (nrow(A) != ncol(A)) stop('powerMethod requires a square matrix')
  if (!all.equal(A, t(A))) stop('powerMethod requires a symmetric matrix')

  d <- nrow(A)
  v <- rnorm(d); v <- v/vector.norm(v)

  for (i in seq_len(maxiter)){
    v_new <- A%*%v; v_new <- v_new/vector.norm(v_new)
    if ((vector.norm(v_new - v)) < eps) break
    v <- v_new
  }
  if (i == maxiter) warning('max iter reached without convergence.')
  return(as.numeric(v_new))
}

