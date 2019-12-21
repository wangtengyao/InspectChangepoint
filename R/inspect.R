inspect <-
function(x, lambda, threshold, schatten = c(1,2), M){
    # inspect is an algorithm for multiple changepoint identification in the 
    # mean structure of a high-dimensional time series
        
    # basic parameters and initialise
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation
    if (missing(lambda)) lambda <- sqrt(log(log(n)*p)/2)
    if (missing(threshold)) threshold <- compute.threshold(n, p)
    if (missing(schatten)) schatten <- 2
    if (missing(M)) M <- 0
    x <- rescale.variance(x)
    
    # generate random time windows of length at least 2
    rnd1 <- sample(0:(n-2), M, replace = TRUE)
    rnd2 <- sample(0:(n-2), M, replace = TRUE)
    window_s <- pmin(rnd1, rnd2)
    window_e <- pmax(rnd1, rnd2) + 2
    
    # recursive function for binary segmentation
    BinSeg <- function(x, s, e, depth, parent.val){
        if (e - s <= 2) return(NULL) # stop when the segment has only one point
        ind <- (window_s >= s) & (window_e <= e) # \mathcal{M}_{s,e}
        max.val <- -1
        cp <- 0

        for (m in c(0,((1:M)[ind]))) {
            if (m == 0) {
                s_m <- s
                e_m <- e
            } else {
                s_m <- window_s[m]
                e_m <- window_e[m]
            }
            obj <- locate.change(x[,(s_m+1):e_m])
            if (obj$cusum > max.val) {
                max.val <- obj$cusum
                cp <- s_m + obj$changepoint
            }
        }
        
        # recurse
        ret <- NULL
        ret$location <- cp
        ret$max.proj.cusum <- max.val #min(parent.val, max.val)
        ret$depth <- depth
        if (ret$max.proj.cusum < threshold) {
            return(NULL) 
        } else {
            return(cbind(BinSeg(x, s, cp, depth+1, ret$max.proj.cusum), 
                         ret, 
                         BinSeg(x, cp, e, depth+1, ret$max.proj.cusum)))
        }
    }
    
    # return all changepoints of x
    ret <- NULL
    ret$x <- x
    ret$changepoints <- BinSeg(x, 0, n, depth = 1, parent.val = .Machine$double.xmax)
    ret$changepoints <- t(matrix(as.numeric(ret$changepoints), nrow = 3))
    colnames(ret$changepoints) = c('location', 'max.proj.cusum', 'depth')
    class(ret) <- 'inspect'
    return(ret)
}
