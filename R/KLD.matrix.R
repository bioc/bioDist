KLD.matrix <- function(x, method=c("locfit", "density"), supp=c(-3,3), 
    subdivisions=1000, diag=FALSE, upper=FALSE){

   x <- as.matrix(x)
   nc <- ncol(x)
   nr <- nrow(x)
   clist <- vector("list", length=nr)
   method = match.arg(method)
   if(method=="locfit")
     {
        for(i in 1:nr)
         clist[[i]] <- locf2func(x[i,])
      }
   else if (method=="density")
      {
         for ( i in 1:nr)
           clist[[i]] <- dens2func(x[i,])
      }
   else 
	stop("method", method, "not supported") 

   rvec<-rep(NA, nr*(nr-1)/2)
   ct <- 1
   for(i in 1:(nr-1))
   for(j in (i+1):nr) {
          if (is.null(supp)) supp <- range(x[c(i,j),])
           rvec[ct] <- KLD(clist[[i]], clist[[j]], supp=supp, 
                                   subdivisions=subdivisions)
           ct <- ct+1
       }
   attributes(rvec) <- list(Size = nr, Labels = row.names(x),
                            Diag = diag, Upper = upper, methods =
                            "KLD", class = "dist")
   rvec

}

