#' Noise standardisation for multivariate time series.
#' @description Each row of the input matrix is normalised by the estimated standard deviation computed through the median absolute deviation of increments.
#' @param x An input matrix of real values.
#' @details This is an auxiliary function used by the \code{InspectChangepoint} package.
#' @return A rescaled matrix of the same size is returned.
#' @examples
#' x <- matrix(rnorm(40),5,8) * (1:5)
#' x.rescaled <- rescale.variance(x)
#' x.rescaled
#' @export

rescale.variance <- function(x){
    p <- dim(x)[1]
    n <- dim(x)[2]
    for (j in 1:p){
        scale <- mad(diff(x[j,]))/sqrt(2)
        x[j,] <- x[j,] / scale
    }
    return(x)
}
