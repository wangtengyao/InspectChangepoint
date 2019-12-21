sparse.svd <-
function(Z, lambda, schatten = c(1, 2), tolerance = 1e-5, max.iter = 10000){
    if (missing(schatten)) schatten <- 2
    if (schatten == 2){
        # with Frobenius norm constraint, the sparse vector is obtained by soft 
        # thresholding
        Mhat <- vector.soft.thresh(Z, lambda)
    } else {
        # with nuclear norm constraint, the sparse vector is obtained by ADMM
        p <- dim(Z)[1]; n <- dim(Z)[2]; gamma <- 1;
        X <- matrix(0,p,n); Y <- matrix(0,p,n); U <- matrix(0,p,n)
        iter <- 0
        while ((iter < max.iter) | (max.iter == 0)) {
            iter <- iter + 1
            X <- PiS(Y - U + gamma * Z)
            Y <- vector.soft.thresh(X + U, lambda * gamma)
            U <- U + (X - Y)
            if (vector.norm(X - Y) < tolerance) break
        }
        Mhat <- X
    }
    vector.proj <- eigen(Mhat%*%t(Mhat), symmetric = TRUE)$vectors[,1]
    return(vector.proj)
}
