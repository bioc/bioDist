"KLdist.matrix" <-
function(x, nbin=10, symmetrize=TRUE, diag=FALSE, upper=FALSE)
{
   x <- t(as.matrix(x))
   nc <- ncol(x)
   nr <- nrow(x)
   clist <- vector("list", length=nr)
    
   appfun <- function(x,y)
    { 
      
      breaks.x <- hist(x,breaks=nbin,freq=FALSE,plot=FALSE)
      breaks.x <- c(breaks.x$breaks[1],breaks.x$breaks[-1][breaks.x$intensities>0])
      
      breaks.y <- hist(y,breaks=nbin,freq=FALSE,plot=FALSE)
      breaks.y <- c(breaks.y$breaks[1],breaks.y$breaks[-1][breaks.y$intensities>0])
      
      temp1 <- table(cut(y,breaks.x))/nc
      ind1 <- temp1>0
      temp1 <- temp1[ind1]
      temp2 <- table(cut(x,breaks.y))/nc
      ind2 <- temp2>0
      temp2 <- temp2[ind2]
      
      dist <- sum(log(temp2)*temp2) - sum(log(table(cut(y,breaks=breaks.y))[ind2>0]/nc)*temp2)
      if(symmetrize)
       {
        dist <- (dist + sum(log(temp1)*temp1)-sum(log(table(cut(x,breaks=breaks.x))[ind1>0]/nc)*temp1))/2
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

}

