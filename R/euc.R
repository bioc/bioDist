setGeneric("euc", function(x, ...) standardGeneric("euc"))

setMethod("euc", signature=signature("matrix"), 
    function(x, diag = FALSE, upper = FALSE)
{
   dist(x, method="euclidean", diag = diag, upper = upper)
} )

setMethod("euc", signature=signature("exprSet"),
    function(x, diag = FALSE, upper = FALSE) euc(x@exprs, diag, upper))

setMethod("euc", signature=signature("ExpressionSet"),
    function(x, diag = FALSE, upper = FALSE) euc(exprs(x), diag, upper))
