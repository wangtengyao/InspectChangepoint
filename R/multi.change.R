#' Generating a high-dimensional time series with multiple changepoints
#' @description The data matrix is generated via X = mu + W, where mu is the mean structure matrix that captures the changepoint locations and sparsity structure, and W is a random noise matrix having independent N(0,sigma^2) entries.
#' @param n Time length of the observation
#' @param p Dimension of the multivariate time series
#' @param ks A vector describing the number of coordinates that undergo a change in each changepoint. If only a scalar is supplied, each changepoint will have the same number of coordinates that undergo a change.
#' @param zs  A vector describing the locations of the changepoints.
#' @param varthetas A vector describing the root mean squared change magnitude in coordinates that undergo a change for each changepoint. If only a scalar is supplied, each changepoint will have the same signal strength value.
#' @param sigma noise level
#' @param overlap A number between 0 and 1. The proportion of overlap in the signal coordinates for successive changepoints.
#' @param shape How the signal strength is distributed across signal coordinates. When shape = 0, all signal coordinates are changed by the same amount; when shape = 1, their signal strength are proportional to 1, sqrt(2), ..., sqrt(k); when shape = 2, they are proportional to 1, 2, ..., k; when shape = 3, they are proportional to 1, 1/sqrt(2), ..., 1/sqrt(k).
#' @return An S3 object of the class 'hdchangeseq' is returned.
#' \itemize{
#'   \item x - The generated data matrix
#'   \item mu - The mean structure of the data matrix
#' }
#'
#' @seealso \code{\link{plot.hdchangeseq}}
#'
#' @examples
#' n <- 2000; p <- 200; ks <- 40;
#' zs <- c(500,1000,1500); varthetas <- c(0.1,0.15,0.2); overlap <- 0.5
#' obj <- multi.change(n, p, ks, zs, varthetas, overlap)
#' plot(obj, noise = TRUE)
#' @export

multi.change <- function(n, p, ks, zs, varthetas, sigma=1, overlap=0, shape=3)
{
    # validating arguments
    nu = length(zs);
    if (length(ks) == 1) {ks = rep(ks, nu)}
    if (length(varthetas) == 1) {varthetas = rep(varthetas, nu)}
    stopifnot(
        nu == length(varthetas),
        nu == length(ks),
        min(zs) >=1,
        max(zs) < n,
        max(ks) <= p
    )

    # generating change vectors
    change = matrix(0, p, nu)
    for (i in 1:nu){
        if (i == 1) {
            start = 0; end = ks[1]-1; ind = (1:ks[1])
        } else if (overlap == -1){
            ind = sample(p,ks[i])
        } else if (overlap == 0){
            start = end + 1; end = start + ks[i] - 1; ind = (start:end)%%p + 1
        } else if (overlap == 1){
            start = 0; end = ks[i] - 1; ind = (start:end)+1;
        } else {
            start = end - round((end - start + 1) * overlap) + 1; end = start + ks[i] - 1; ind = (start:end)%%p + 1
        }

        if (shape == 0) {
            theta = rep(1,ks[i]); theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        } else if (shape == 1) {
            theta = ks[i]:1; theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        } else if (shape == 2){
            theta = (ks[i]:1)^2; theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        } else if (shape == 3){
            theta = (1:ks[i])^(-1/2); theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        }
        change[ind, i] = theta
    }

    # generating mean structure
    mu <- matrix(0, p, n)
    for (i in 1:nu) {mu[,(zs[i]+1):n] = mu[,(zs[i]+1):n] + change[,i]}

    # generating data matrix
    x <- mu + matrix(rnorm(p*n),p,n)*sigma
    ret <- NULL
    ret$x = x; ret$mu = mu;
    class(ret) <- 'hdchangeseq'
    return(ret);
}
