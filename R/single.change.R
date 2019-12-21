single.change <-
function(n, p, k, z, vartheta, sigma = 1, shape = 3, noise = 0, corr = 0){
    mu <- matrix(0, p, n)
    if (shape == 0) {
        theta = rep(vartheta, k)
    } else if (shape == 1) {
        theta = vartheta * sqrt((1:k)/((k+1)/2))
    } else if (shape == 2){
        theta = vartheta * (1:k) / sqrt((k+1)*(2*k+1)/6)      
    } else if (shape == 3){
        theta = (1:k)^(-1/2); theta = theta/vector.norm(theta)*sqrt(k)*vartheta
    }
    if (noise != -1) {
        mu[1:k, (z+1):n] <- theta
    } else { # non-simultaneous changepoints
        for (i in 1:k){b = sample(-corr:corr,1); mu[i,(z+b):n] <- theta[i]}
    }
    
    if (noise <= 0) { # independent normal noise
        W = matrix(rnorm(p*n), p, n)* sigma
    } else if (noise == 0.5){ # local cross-sectional correlation
        W = t(
            mvrnorm(n, mu = rep(0,p), Sigma = outer(1:p,1:p,function(a,b){corr^(abs(a-b))}))) * sigma
    } else if (noise == 1){ # global cross-sectional correlation
        W = t(mvrnorm(n, mu = rep(0,p), Sigma = matrix(corr,p,p)+diag(p)*(1-corr))) * sigma
    } else if (noise == 2){ # temporal AR noise
        W = matrix(0,p,n); W[,1] = rnorm(p)
        for (j in 1:p){
            for (t in 2:n) W[j,t] = W[j,t-1]*sqrt(corr) + rnorm(1)*sqrt(1-corr)
        }
        W = W * sigma
    } else if (noise == 3){ # uniform noise
        W = matrix(runif(p*n, -sqrt(3),sqrt(3)), p, n)*sigma
    } else if (noise == 4){ # centred exponential noise
        W = matrix(rexp(p*n)-1, p, n)*sigma
    } 

    x <- mu +  W 
    ret <- NULL
    ret$x <- x; ret$mu <- mu; 
    class(ret) <- 'hdchangeseq'
    return(ret)
}
