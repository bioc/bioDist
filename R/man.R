man <- function(x, diag = FALSE, upper = FALSE)
{
  dist(x, method="manhattan", diag=diag, upper=upper)
}

