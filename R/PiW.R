#' Projection onto the standard simplex
#' @description The input vector is projected onto the standard simplex, i.e. the set of vectors of the same length as the input vector with non-negative entries that sum to 1.
#' @param v Input vector
#' @details This is an auxiliary function used by the \code{InspectChangepoint} package.
#' @return A vector in the standard simplex that is closest to the input vector is returned.
#' @references Chen, Y. and Ye, X. (2011) Projection onto a simplex. arXiv preprint, arxiv:1101.6081.
#' @examples
#' v <- rnorm(10)
#' PiW(v)
#' @export

PiW <- function(v){
    # projection of a vector to the standard simplex (non-negative entries
    # summing to 1)
    ord = order(v, decreasing = TRUE);
    v = v[ord]
    s = 0
    for (i in 1:length(v)){
        s = s + v[i]
        if (i == length(v)) {delta = (s-1)/i; break}
        d = s - v[i+1]*i
        if (d >= 1) {delta = (s - 1)/i; break}
    }
    v[1:i] = v[1:i] - delta; v[-(1:i)] = 0;
    w = rep(0, length(v))
    w[ord] = v
    w
}
