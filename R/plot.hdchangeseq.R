#' Plot function for 'hdchangeseq' class
#' @description Visualising the high-dimensional time series in an 'hdchangeseq' class object. The data matrix or its mean structure is visualised using a grid of coloured rectangles with colours corresponding to the value contained in corresponding coordinates. A heat-spectrum (red to white for values low to high) is used to convert values to colours.
#'
#' @param x An object of 'hdchangeseq' class
#' @param noise If noise == TRUE, the data matrix is plotted, otherwise, only the mean structure is plotted.
#' @param shuffle Whether to shuffle the rows of the plotted matrix.
#' @param ... Other graphical parameters are not used.
#'
#' @examples
#' n <- 2000; p <- 200; ks <- 40; zs <- c(500,1000,1500)
#' varthetas <- c(0.1,0.15,0.2); overlap <- 0.5
#' obj <- multi.change(n, p, ks, zs, varthetas, overlap)
#' plot(obj, noise = TRUE)
#' @export

plot.hdchangeseq <- function(x, noise = TRUE, shuffle = FALSE, ...){
    obj <- x
    p <- dim(obj$x)[1]
    if (shuffle) {
      ind = sample(p); obj$x = obj$x[ind,]; obj$mu = obj$mu[ind,]
    }
    if (noise == TRUE) {
      image(t(obj$x[p:1,]), axes = FALSE, frame.plot = FALSE)
    } else {
      image(t(obj$mu[p:1,]), axes = FALSE, frame.plot = FALSE)
    }
}
