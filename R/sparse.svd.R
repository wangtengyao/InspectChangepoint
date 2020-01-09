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
#' @param tolerance Tolerance criterion for convergence of the ADMM algorithm. Not used when shatten=2.
#' @param max.iter Maximum number of iteration in the ADMM algorithm. Not used when shatten=2.
#'
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

    if (nrow(Mhat) < ncol(Mhat)){
        vector.proj <- power.method(Mhat%*%t(Mhat), 1e-5)
    } else {
        tmp <- Mhat %*% power.method(t(Mhat)%*%Mhat, 1e-5)
        vector.proj <- tmp/vector.norm(tmp)
    }

    return(vector.proj)
}

