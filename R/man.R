setGeneric("man", function(x, ...) standardGeneric("man"))

setMethod("man", signature=signature("matrix"),
    function(x, diag = FALSE, upper = FALSE)
{
  dist(x, method="manhattan", diag=diag, upper=upper)
} )

setMethod("man", signature=signature("exprSet"),
    function(x, diag=FALSE, upper=FALSE) man(x@exprs, diag, upper))

setMethod("man", signature=signature("ExpressionSet"),
    function(x, diag=FALSE, upper=FALSE) man(exprs(x), diag, upper))

