"man" <-
function(x, diag = FALSE, upper = FALSE)
{
  dist(t(x), method="manhattan", diag=diag, upper=upper)
}

