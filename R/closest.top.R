"closest.top" <-
function(x, dist.mat, top)
{
  dist <- as.matrix(dist.mat)
  vector <- dist[x,colnames(dist) != x]
  return(names(vector)[order(vector)[1:top]])
}

