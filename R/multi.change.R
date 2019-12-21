multi.change <-
function(n, p, ks, zs, varthetas, sigma = 1, overlap = 0, shape = 3)
    # generate a multivariate time series panel data with n time points and p dimensional
    # zs varthetas and ks are equal length vectors describing location, 
    # magnitude and sparsity of changepoints
    # overlap describes how the changepoint coordinates are chosen, 
    # 0 means no overlap, 1 means maximal overlap
    # 0.5 means half overlap, -1 means random
{
    # validating arguments
    nu = length(zs); if (length(ks) == 1) {ks = rep(ks, nu)}; if (length(varthetas) == 1) {varthetas = rep(varthetas, nu)}
    stopifnot(nu == length(varthetas), nu == length(ks), min(zs) >=1, max(zs) < n, max(ks) <= p)
    
    # generating change vectors
    change = matrix(0, p, nu)
    for (i in 1:nu){
        if (i == 1) {
            start = 0; end = ks[1]-1; ind = (1:ks[1])
        } else if (overlap == -1){
            ind = sample(p,ks[i])
        } else if (overlap == 0){
            start = end + 1; end = start + ks[i] - 1; ind = (start:end)%%p + 1
        } else if (overlap == 1){
            start = 0; end = ks[i] - 1; ind = (start:end)+1; 
        } else {
            start = end - round((end - start + 1) * overlap) + 1; end = start + ks[i] - 1; ind = (start:end)%%p + 1
        }
        
        if (shape == 0) {
            theta = rep(1,ks[i]); theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        } else if (shape == 1) {
            theta = ks[i]:1; theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        } else if (shape == 2){
            theta = (ks[i]:1)^2; theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        } else if (shape == 3){
            theta = (1:ks[i])^(-1/2); theta = theta/vector.norm(theta)*sqrt(ks[i])*varthetas[i]
        }
        change[ind, i] = theta
    }
    
    # generating mean structure
    mu <- matrix(0, p, n)
    for (i in 1:nu) {mu[,(zs[i]+1):n] = mu[,(zs[i]+1):n] + change[,i]}

    # generating data matrix
    x <- mu + matrix(rnorm(p*n),p,n)*sigma
    ret <- NULL
    ret$x = x; ret$mu = mu; 
    class(ret) <- 'hdchangeseq'
    return(ret);
}
