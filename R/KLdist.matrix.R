setGeneric("KLdist.matrix", function(x, ...) standardGeneric("KLdist.matrix"))

setMethod("KLdist.matrix", signature=signature("matrix"), 
          function(x, nbin=10, symmetrize=FALSE, diag=FALSE, upper=FALSE)
{
   x <- as.matrix(x)
   nc <- ncol(x)
   nr <- nrow(x)
   clist <- vector("list", length=nr)
   me <- .Machine$double.eps
 
   ##note: we combine x and y before binning, to make sure we span
   ##   the range of the data, and we add machine epsilon to 
   ##   protect against +/- Inf; this could use some work.
   appfun <- function(x,y)
    { 
      
      breaks.x <- hist(c(x,y) ,breaks=nbin,plot=FALSE)$breaks
      
      temp1 <- table(cut(y,breaks.x))/nc
      temp1 <- temp1+me
      temp2 <- table(cut(x,breaks.x))/nc
      temp2 <- temp2 + me
      
      dist <- sum(log(temp2/temp1)*temp2) 
      if(symmetrize)
       {
        dist <- (dist + sum(log(temp1/temp2)*temp1))/2
         }
      return(dist)   
    }
    
   rvec<-rep(NA, nr*(nr-1)/2)
   ct <- 1
   for(i in 1:(nr-1))
       for(j in (i+1):nr) {
           rvec[ct] <- appfun(x[i,], x[j,])
           ct <- ct+1
       }
   attributes(rvec) <- list(Size = nr, Labels = row.names(x),
                            Diag = diag, Upper = upper, methods =
                            "KLdist", class = "dist")
   rvec

} )

setMethod("KLdist.matrix", signature=signature("exprSet"),
function(x, nbin=10, symmetrize=FALSE, diag=FALSE, upper=FALSE, 
         sample=TRUE) {
    .Deprecated(msg=EXPRSET_MSG)
    if( sample ) ep = t(exprs(x)) else ep = exprs(x)
    KLdist.matrix(ep, nbin, symmetrize, diag, upper)
})

setMethod("KLdist.matrix", signature=signature("ExpressionSet"),
function(x, nbin=10, symmetrize=FALSE, diag=FALSE, 
         upper=FALSE, sample=TRUE)  {
    if( sample ) ep = t(exprs(x)) else ep = exprs(x)
    KLdist.matrix(ep, nbin, symmetrize, diag, upper)
})




## tentative "list" method for unequal sized samples (added by
## Deepayan Sarkar)


setMethod("KLdist.matrix",
          signature=signature("list"), 
          function(x, nbin=10, symmetrize=FALSE,
                   diag=FALSE, upper=FALSE)
      {
          n <- length(x)
          clist <- vector("list", length=n)
          me <- .Machine$double.eps
          
          ##note: we combine x and y before binning, to make sure we span
          ##   the range of the data, and we add machine epsilon to 
          ##   protect against +/- Inf; this could use some work.
          appfun <- function(x,y)
          { 
              
              breaks.x <- hist(c(x,y) ,breaks=nbin,plot=FALSE)$breaks
              
              temp1 <- table(cut(y,breaks.x))/nc
              temp1 <- temp1+me
              temp2 <- table(cut(x,breaks.x))/nc
              temp2 <- temp2 + me
              
              dist <- sum(log(temp2/temp1)*temp2) 
              if(symmetrize)
              {
                  dist <- (dist + sum(log(temp1/temp2)*temp1))/2
              }
              return(dist)   
          }
          ans <- matrix(NA, n, n)
          for(i in seq_len(n))
              for(j in seq_len(n))
              {
                  ans[i, j] <- appfun(x[[i]], x[[j]])
              }
          ans
      } )






