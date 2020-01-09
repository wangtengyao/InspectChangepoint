#' Soft thresholding a vector
#' @param x a vector of real numbers
#' @param lambda soft thresholding value
#' @return a vector of the same length
#' @description entries of v are moved towards 0 by the amount lambda until they hit 0.
#' @export

vector.soft.thresh <- function(x, lambda){
  sign(x)*pmax(0,(abs(x) - lambda))
}
