setGeneric("cor.dist", function(x, ...) standardGeneric("cor.dist"))

setMethod("cor.dist", signature=signature("matrix"), 
    function(x, abs=TRUE,diag=FALSE, upper=FALSE)
{
  nr <- nrow(x)
  rvec <- cor(t(x))
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
                              "cor", class = "dist")
   rvec
} )

setMethod("cor.dist", signature=signature("eSet"),
    function(x, abs=TRUE,diag=FALSE, upper=FALSE) {
        if( sample ) ep = t(exprs(x)) else ep = exprs(x)
        cor.dist(ep, abs, diag, upper)
    })



