setGeneric("tau.dist", function(x, ...) standardGeneric("tau.dist"))

setMethod("tau.dist", signature=signature("matrix"), 
    function(x, abs=TRUE,diag=FALSE, upper=FALSE)
{
  nr <- nrow(x)
  rvec <- cor(t(x), method="kendall")
  if(abs)
   rvec <- 1-abs(rvec)
  else
    rvec <- 1-rvec
  if(upper)
   rvec <- rvec[upper.tri(rvec,diag=diag)]
  else
     rvec <- rvec[lower.tri(rvec,diag=diag)]
  attributes(rvec) <- list(Size = nr, Labels = rownames(x),
                              Diag = diag, Upper = upper, methods =
                              "kendall", class = "dist")
   rvec
} )

setMethod("tau.dist", signature=signature("exprSet"),
    function(x, abs=TRUE,diag=FALSE, upper=FALSE) {
        .Deprecated(msg=EXPRSET_MSG)
        tau.dist(x@exprs, abs, diag, upper)
        })

setMethod("tau.dist", signature=signature("ExpressionSet"),
    function(x, abs=TRUE,diag=FALSE, upper=FALSE) 
        tau.dist(exprs(x), abs, diag, upper))
