"cor.dist" <-
function(x, abs=TRUE,diag=FALSE, upper=FALSE)
{
  nc <- ncol(x)
  rvec <- cor(x)
  if(abs)
   rvec <- 1-abs(rvec)
  else
    rvec <- 1-rvec
  if(upper)
   rvec <- rvec[upper.tri(rvec,diag=diag)]
  else
     rvec <- rvec[lower.tri(rvec,diag=diag)]
  attributes(rvec) <- list(Size = nc, Labels = colnames(x),
                              Diag = diag, Upper = upper, methods =
                              "cor", class = "dist")
   rvec
}

