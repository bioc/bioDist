"euc" <-
function(x, diag = FALSE, upper = FALSE)
{
   dist(t(x), method="euclidean", diag = diag, upper = upper)
}

