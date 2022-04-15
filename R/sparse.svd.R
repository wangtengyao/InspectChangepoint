#' Computing the sparse leading left singular vector of a matrix
#'
#' @description Estimating the sparse left leading singular vector by first computing a maximiser Mhat of the convex problem
#' \deqn{<Z, M> - \lambda |M|_1}
#' subject to the Schatten norm constraint |M|_schatten <= 1 using alternating direction method of multipliers (ADMM). Then the leading left singular vector of Mhat is returned.
#'
#' @details In case of schatten = 2, a closed-form solution for Mhat using matrix soft thresholding is possible. We use the closed-form solution instead of the ADMM algorithm to speed up the computation.
#'
#' @param Z Input matrix whose left leading singular vector is to be estimated.
#' @param lambda Regularisation parameter
#' @param schatten Schatten norm constraint to be used. Default uses Schatten-2-norm, i.e. the Frobenius norm. Also possible to use Schatten-1-norm, the nuclear norm.
#' @return A vector that has the same length as nrow(Z) is returned.
#' @examples
#' Z <- matrix(rnorm(20),4,5)
#' lambda <- 0.5
#' sparse.svd(Z, lambda)
#' @export

sparse.svd <- function(Z, lambda, schatten=c(1, 2), tolerance=1e-5, max.iter=10000){
  if (missing(schatten)) schatten <- 2
  if (schatten == 2){
    # with Frobenius norm constraint, the sparse vector is obtained by soft
    # thresholding
    Mhat <- vector.soft.thresh(Z, lambda)
  } else {
    # with nuclear norm constraint, the sparse vector is obtained by ADMM
    p <- dim(Z)[1]; n <- dim(Z)[2]; gamma <- 1;
    X <- matrix(0,p,n); Y <- matrix(0,p,n); U <- matrix(0,p,n)
    iter <- 0
    while ((iter < max.iter) | (max.iter == 0)) {
      iter <- iter + 1
      X <- PiS(Y - U + gamma * Z)
      Y <- vector.soft.thresh(X + U, lambda * gamma)
      U <- U + (X - Y)
      if (vector.norm(X - Y) < tolerance) break
    }
    Mhat <- X
  }

  # compute the leading left singular vector of Mhat
  if (sum(Mhat^2)!=0){
    # compute the leading left singular vector
    if (require(RSpectra, attach.required=FALSE) && min(dim(Mhat)) >= 3){
      vector.proj <- RSpectra::svds(Mhat, 1)$u[,1]
    } else {
      vector.proj <- svd(Mhat)$u[,1]
    }
  } else {
    # if the thresholded matrix is zero, return a random vector
    vector.proj <- rnorm(p)
    vector.proj <- vector.proj / vector.norm(vector.proj)
  }
  return(vector.proj)
}

#' Computing the sparse leading left singular vector of a matrix with missing entries
#' @param Z Input matrix whose left leading singular vector is to be estimated.
#' @param lambda Regularisation parameter
#' @param max_iter maximum iteration
#' @param tol tolerance level for convergence
sparse.svd.missing <- function(Z, lambda, max_iter=1000, tol=1e-10){

  if (sum(abs(Z)) == 0) return(random.UnitVector(nrow(Z)))

  if (require(RSpectra, attach.required=FALSE) && min(dim(Z)) >= 3){
    vhat <- RSpectra::svds(Z, 1)$u[,1]
  } else {
    vhat <- svd(Z)$u[,1]
  }

  for (iter in 1:max_iter){
    vhat_old <- vhat
    what <- vector.normalise(t(Z) %*% vhat)
    tmp <- Z %*% what
    lambda_tmp <- min(lambda, max(abs(tmp)) - 1e-10)
    vhat <- vector.normalise(vector.soft.thresh(Z %*% what, lambda_tmp))
    if (vector.norm(vhat_old - vhat) < tol) break
  }

  return(as.vector(vhat))
}
