PiS <-
function(M){
    # projection of a matrix (in Frobenius norm) to the ball of nuclear norm 1
    tmp = svd(M); 
    tmp$u%*%diag(PiW(tmp$d))%*%t(tmp$v)
}
