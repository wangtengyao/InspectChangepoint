#' Computing threshold used in \code{inspect}
#' @description The threshold level to be used in \code{inspect} is computed via Monte Carlo simulation of multivariate time series that do not contain any changepoints.
#' @param n Time length of the observation.
#' @param p Dimension of the multivariate time series.
#' @param nrep Number of Monte Carlo repetition to be used.
#' @param show_progress whether to show the progress of Monte Carlo simulation
#' @return A numeric value indicating the threshold level that should be used based on the Monte Carlo simulation.
#' @examples
#' compute.threshold(n=200, p=50)
#' @export

compute.threshold <- function(n, p, nrep=100, show_progress=TRUE){
    if (show_progress) cat('Calculating threshold... ')
    cusum.stats = rep(0,nrep)
    for (i in 1:nrep) {
        x <- single.change(n, p, 1, n - 1, 0)$x
        cusum.stats[i] = locate.change(x)$cusum
        if (show_progress) printPercentage(i, nrep)
    }
    thresh <- as.numeric(max(cusum.stats))
    cat('. Threshold =', thresh, '\n')
    return(thresh)
}
