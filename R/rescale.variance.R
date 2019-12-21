rescale.variance <-
function(x){
    # standardise each row of the matrix x by median absolute deviation
    p <- dim(x)[1]
    n <- dim(x)[2]
    for (j in 1:p){
        scale <- mad(diff(x[j,]))/sqrt(2)
        x[j,] <- x[j,] / scale
    }
    return(x)
}
