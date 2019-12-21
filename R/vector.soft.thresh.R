vector.soft.thresh <-
function(x, lambda){
    # soft thresholding, moving entries of x towards 0 by amount lambda,
    # until hits 0
    sign(x)*pmax(0,(abs(x)-lambda))
}
