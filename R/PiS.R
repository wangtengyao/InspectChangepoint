#' Matrix projection onto the nuclear norm unit sphere
#' @description Projection (with respect to the inner product defined by the Frobenius norm) of a matrix onto the unit sphere defined by the nuclear norm.
#' @details This is an auxiliary function used by the \code{InspectChangepoint} package. The projection is achieved by first performing a singular value decomposition, then projecting the vector of singular values onto the standard simplex, and finally using singular value decomposition in reverse to build the projected matrix.
#' @param M Input matrix
#' @return A matrix of the same dimension as the input is returned.
#' @examples
#' M <- matrix(rnorm(20),4,5)
#' PiS(M)
#' @export

PiS <- function(M){
    # projection of a matrix (in Frobenius norm) to the ball of nuclear norm 1
    tmp = svd(M);
    tmp$u%*%diag(PiW(tmp$d))%*%t(tmp$v)
}
