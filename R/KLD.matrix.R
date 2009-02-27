setGeneric("KLD.matrix", function(x, ...) standardGeneric("KLD.matrix"))

setMethod("KLD.matrix", signature=signature("matrix"), 
    function(x, method=c("locfit", "density"), supp=c(-3,3), 
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
          if (is.null(supp)) supp <- extendrange(x[c(i,j),],r=range(x,na.rm=T),f=0.01)
           rvec[ct] <- KLD(clist[[i]], clist[[j]], supp=supp,subdivisions=subdivisions)
           ct <- ct+1
       }
   attributes(rvec) <- list(Size = nr, Labels = row.names(x),
                            Diag = diag, Upper = upper, methods =
                            "KLD", class = "dist")
   rvec

})

setMethod("KLD.matrix", signature=signature("ExpressionSet"),
    function(x, method=c("locfit", "density"), supp=c(-3,3), 
        subdivisions=1000, diag=FALSE, upper=FALSE, sample=TRUE) {
        if( sample ) ep = t(exprs(x)) else ep = exprs(x)
        KLD.matrix(ep, method, supp, subdivisions, diag, upper)
        })


## tentative "list" method for unequal sized samples (added by
## Deepayan Sarkar)

setMethod("KLD.matrix", signature=signature("list"), 
          function(x, method = c("locfit", "density"),
                   supp=c(-3,3), 
                   subdivisions=1000,
                   diag=FALSE, upper=FALSE)
      {  
          method <- match.arg(method)
          dfun <- switch(method,
                         locfit = locf2func,
                         density = dens2func)
          n <- length(x)
          if (n < 1) return()
          clist <- vector("list", length=n)
          #for (i in seq_len(n)) clist[[i]] <- dfun(x[[i]])
	  for (i in seq_len(n)) clist[[i]] <- dens2func(x[[i]])
 
	    
	  rvec<-rep(NA, n*(n-1)/2)
	      ct <- 1
	      for(i in 1:(n-1)){
	      for(j in (i+1):n) {
                      if (is.null(supp)) supp <- range(c(x[[i]],x[[j]]))
		      rvec[ct] <- KLD(clist[[i]], clist[[j]], supp=supp,subdivisions=subdivisions)
		      ct <- ct+1
		  }
               } 
              attributes(rvec) <- list(Size = n, Labels = row.names(x),
					Diag = diag, Upper = upper, methods =
					"KLD", class = "dist")
	     
              return(rvec)

#           ans <- matrix(NA, n, n)
#           for(i in seq_len(n))
#               for(j in seq_len(n))
#               {
#                   if (is.null(supp))
#                       supp <- range(x[[i]], x[[j]], finite = TRUE)
#                   ans[i, j] <-
#                       if (i == j) 0
#                       else KLD(clist[[i]], clist[[j]],
#                                supp=supp, 
#                                subdivisions=subdivisions)
#               }
#           ## if (symmetrize) ans <- t(ans) + ans 
#           if (!is.null(names(x)))
#               rownames(ans) <- colnames(ans) <- names(x)
#           ans


      })

