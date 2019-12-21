plot.hdchangeseq <-
function(x, noise = TRUE, shuffle = FALSE, ...){
    obj <- x
    p = dim(obj$x)[1]
    if (shuffle) {ind = sample(p); obj$x = obj$x[ind,]; obj$mu = obj$mu[ind,]}
    if (noise == TRUE) image(t(obj$x[p:1,]), axes = FALSE, frame.plot = FALSE)
    else image(t(obj$mu[p:1,]), axes = FALSE, frame.plot = FALSE)
}
