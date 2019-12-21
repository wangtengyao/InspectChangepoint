locate.change <-
function(x, lambda, schatten = 2, sample.splitting = FALSE, standardize.series = FALSE) 
# find a single changepoint in multivariate time series x
{
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation
    if (missing(lambda)) lambda <- sqrt(log(log(n)*p)/2)
    if (standardize.series) x <- rescale.variance(x)
    if (sample.splitting){
        x1 <- x[,seq(1,n,by=2)]
        x2 <- x[,seq(2,n,by=2)]
    } else {
        x1 <- x
        x2 <- x
    }
    
    # construct cusum matrix of x
    cusum.matrix1 <- cusum.transform(x1)
    if (sample.splitting) cusum.matrix2 <- cusum.transform(x2) else cusum.matrix2 <- cusum.matrix1
    
    # estimate changepoint
    if (lambda >= max(abs(cusum.matrix1))) lambda <- max(abs(cusum.matrix1)) - 1e-10
    
    vector.proj <- sparse.svd(cusum.matrix1, lambda, schatten);
    cusum.proj <- t(cusum.matrix2)%*%vector.proj
    
    ret <- NULL
    ret$changepoint <- which.max(abs(cusum.proj))
    if (sample.splitting) ret$changepoint <- ret$changepoint * 2
    ret$cusum <- max(abs(cusum.proj))
    return(ret)
}
