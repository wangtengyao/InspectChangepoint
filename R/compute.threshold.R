compute.threshold <-
function(n,p,nrep = 100){
    cusum.stats = rep(0,nrep)
    for (i in 1:nrep) {
        cusum.stats[i] = locate.change(single.change(n,p,1,n-1,0)$x)$cusum
    }
    as.numeric(max(cusum.stats))
}
