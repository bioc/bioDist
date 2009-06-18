setGeneric("euc", function(x, ...) standardGeneric("euc"))

setMethod("euc", signature=signature("matrix"), 
    function(x, diag = FALSE, upper = FALSE)
{
   dist(x, method="euclidean", diag = diag, upper = upper)
} )

setMethod("euc", signature=signature("eSet"),
    function(x, diag = FALSE, upper = FALSE,sample =TRUE) {
    if(sample) ep <- t(exprs(x)) else ep <- exprs(X)
    euc(ep, diag, upper)
})
