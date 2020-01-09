#' Generating high-dimensional time series with exactly one change in the mean
#' structure
#' @description The data matrix is generated via X = mu + W, where mu is the mean structure matrix that captures the changepoint location and sparsity structure, and W is a random noise matrix.
#'
#' @param n Time length of the observation
#' @param p Dimension of the multivariate time series
#' @param k Number of coordinates that undergo a change
#' @param z Changepoint location, a number between 1 and n-1.
#' @param vartheta The root mean squared change magnitude in coordinates that undergo a change
#' @param sigma noise level, see \code{noise} for more details.
#' @param shape How the signal strength is distributed across signal coordinates. When shape = 0, all signal coordinates are changed by the same amount; when shape = 1, their signal strength are proportional to 1, sqrt(2), ..., sqrt(k); when shape = 2, they are proportional to 1, 2, ..., k; when shape = 3, they are proportional to 1, 1/sqrt(2), ..., 1/sqrt(k).
#' @param noise Noise structure of the multivarite time series. For noise = 0, 0.5, 1, columns of W have independent multivariate normal distribution with covariance matrix Sigma. When noise = 0, Sigma = sigma^2 * I_p; when noise = 0.5, noise has local dependence structure given by Sigma_{i,j} = sigma*corr^|i-j|; when noise = 1, noise has global dependence structure given by matrix(corr,p,p)+diag(p)*(1-corr))) * sigma. When noise = 2, rows of the W are independent and each having an AR(1) structure given by W_{j,t} = W_{j,t-1} * sqrt(corr) + rnorm(sd = sigma) * sqrt(1-corr). For noise = 3, 4, entries of W have i.i.d. uniform distribution and exponential distribution respectively, each centred and rescaled to have zero mean and variance sigma^2.
#' @param corr Used to specify correlation structure in the noise. See \code{noise} for more details.
#'
#' @return An S3 object of the class 'hdchangeseq' is returned.
#' \itemize{
#'   \item x - The generated data matrix
#'   \item mu - The mean structure of the data matrix
#' }
#'
#' @seealso \code{\link{plot.hdchangeseq}}
#'
#' @examples
#' n <- 2000; p <- 100; k <- 10; z <- 800; vartheta <- 1; sigma <- 1
#' shape <- 3; noise <- 0; corr <- 0
#' obj <- single.change(n,p,k,z,vartheta,sigma, shape, noise, corr)
#' plot(obj, noise = TRUE)
#' @importFrom MASS mvrnorm
#' @export

single.change <- function(n, p, k, z, vartheta, sigma = 1, shape = 3, noise = 0,
                          corr = 0){
    mu <- matrix(0, p, n)
    if (shape == 0) {
        theta = rep(vartheta, k)
    } else if (shape == 1) {
        theta = vartheta * sqrt((1:k)/((k+1)/2))
    } else if (shape == 2){
        theta = vartheta * (1:k) / sqrt((k+1)*(2*k+1)/6)
    } else if (shape == 3){
        theta = (1:k)^(-1/2); theta = theta/vector.norm(theta)*sqrt(k)*vartheta
    }
    if (noise != -1) {
        mu[1:k, (z+1):n] <- theta
    } else { # non-simultaneous changepoints
        for (i in 1:k){b = sample(-corr:corr,1); mu[i,(z+b):n] <- theta[i]}
    }

    if (noise <= 0) { # independent normal noise
        W = matrix(rnorm(p*n), p, n)* sigma
    } else if (noise == 0.5){ # local cross-sectional correlation
        W = t(
            mvrnorm(n, mu = rep(0,p), Sigma = outer(1:p,1:p,function(a,b){corr^(abs(a-b))}))) * sigma
    } else if (noise == 1){ # global cross-sectional correlation
        W = t(mvrnorm(n, mu = rep(0,p), Sigma = matrix(corr,p,p)+diag(p)*(1-corr))) * sigma
    } else if (noise == 2){ # temporal AR noise
        W = matrix(0,p,n); W[,1] = rnorm(p)
        for (j in 1:p){
            for (t in 2:n) W[j,t] = W[j,t-1]*sqrt(corr) + rnorm(1)*sqrt(1-corr)
        }
        W = W * sigma
    } else if (noise == 3){ # uniform noise
        W = matrix(runif(p*n, -sqrt(3),sqrt(3)), p, n)*sigma
    } else if (noise == 4){ # centred exponential noise
        W = matrix(rexp(p*n)-1, p, n)*sigma
    }

    x <- mu +  W
    ret <- NULL
    ret$x <- x; ret$mu <- mu;
    class(ret) <- 'hdchangeseq'
    return(ret)
}
