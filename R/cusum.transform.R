cusum.transform <-
function(x){
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation
    
    leftsums <- t(apply(x,1,cumsum))
    rightsums <- leftsums[,n]-leftsums
    leftsums <- t(leftsums)
    rightsums <- t(rightsums)
    t <- 1:(n-1)
    
    # constructing CUSUM matrix
    return(t((rightsums[t,]/(n-t) - leftsums[t,]/t)*sqrt(t*(n-t)/n)))
}
