PiW <-
function(v){
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
