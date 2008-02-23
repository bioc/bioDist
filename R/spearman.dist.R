setGeneric("spearman.dist", function(x, ...) standardGeneric("spearman.dist"))

setMethod("spearman.dist", signature=signature("matrix"), 
    function(x, abs=TRUE,diag=FALSE, upper=FALSE)
{
  nr <- nrow(x)
  rvec <- cor(t(x), method="spearman")
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
                              "spearman", class = "dist")
   rvec
} )

setMethod("spearman.dist", signature=signature("ExpressionSet"),
    function(x, abs=TRUE,diag=FALSE, upper=FALSE, sample=TRUE) {
        if( sample ) ep = t(exprs(x)) else ep = exprs(x)
        spearman.dist(ep, abs, diag, upper)})
