plot.inspect <-
function(x, ...){
    obj <- x
    p = dim(obj$x)[1]
    n = dim(obj$x)[2]
    image(1:n, 1:p, t(obj$x), xlab = 'time', ylab = 'coordinate', axes = TRUE, frame.plot = TRUE)
    result = as.matrix(obj$changepoints)
    abline(v = result[,1], col = 'blue')
}
