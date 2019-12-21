vector.norm <-
function(v, q = 2){
    # calculate the L^q norm of a vector or a matrix
    if (q > 0) (sum(abs(v)^q))^(1/q)
    else if (q == 0) sum(v!=0)
    else if (q == -1) max(abs(v))
    else NaN
}
