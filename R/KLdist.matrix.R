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

   interpfunc <- function(x,y,...) {
	f <- function(w) approx(x,y, w, yleft=0, yright=0)$y
	class(f) <- "dfun"
	f
   }
   #### to be optimized  :ngopalak
   appfun <- function(x,y)
    { 
	h1 <- dpih(x)
 	bins <- seq(min(x)-0.1, max(x)+0.1+h1, by=h1)
        temp1Hist <- hist(x,breaks=bins,plot=FALSE)
        temp1 <- temp1Hist$count/(nc*h1)

	h2 <- dpih(y)
 	bins <- seq(min(y)-0.1, max(y)+0.1+h2, by=h2)     
 	temp2Hist <- hist(y,breaks=bins,plot=FALSE)
 	temp2 <- temp2Hist$count/(nc*h2)
 	temp2 <- temp2
        
        f <- interpfunc(temp1Hist$mids,temp1)
        g<- interpfunc(temp2Hist$mids,temp2)
        step <- min(c(h1,h2))/2
        supp <- range(c(x,y))
        p<- seq(from= supp[1], to =supp[2], by= step)        
        dist<-sum(log((f(p)+me)/(g(p)+me))*f(p))*step
        if(symmetrize)
        {
          dist <- (dist +  sum(log((g(p)+me)/(f(p)+me))*g(p))*step)/2

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


setMethod("KLdist.matrix", signature=signature("ExpressionSet"),
function(x, nbin=10, symmetrize=FALSE, diag=FALSE, 
         upper=FALSE, sample=TRUE)  {
    if( sample ) ep = t(exprs(x)) else ep = exprs(x)
    KLdist.matrix(ep, nbin, symmetrize, diag, upper)
})


setMethod("KLdist.matrix",
          signature=signature("list"), 
          function(x, 
                   discretize = TRUE, nbin = 10,
                   symmetrize = FALSE,
                   diag = FALSE, upper=FALSE)
      {
          n <- length(x)
          clist <- vector("list", length=n)
          me <- .Machine$double.eps
          interpfunc <- function(x,y,...) {
            f <- function(w) approx(x,y, w, yleft=0, yright=0)$y
            class(f) <- "dfun"
            f
          }
          ##note: we combine x and y before binning, to make sure we span
          ##   the range of the data, and we add machine epsilon to 
          ##   protect against +/- Inf; this could use some work.
          distfun <- function(x, y)
          {   
              ## not clear what should be done if exactly one of x and y is a factor
              if (discretize && !is.factor(x))
              {  
                  ## to be optimized , use table instead of hist
	          h1 <- dpih(x)
		  nc1<- length(x)
 		  bins <- seq(min(x)-0.1, max(x)+0.1+h1, by=h1)
        	  temp1Hist <- hist(x,breaks=bins,plot=FALSE)
                  temp1 <- temp1Hist$count/(nc1*h1)
		 
                  h2 <- dpih(y)
		  nc2 <-length(y)
 	          bins <- seq(min(y)-0.1, max(y)+0.1+h2, by=h2)     
 	          temp2Hist <- hist(y,breaks=bins,plot=FALSE)
 	          temp2 <- temp2Hist$count/(nc2*h2)
 	          temp2 <- temp2
                  
         	  f <- interpfunc(temp1Hist$mids,temp1)
                  g<- interpfunc(temp2Hist$mids,temp2)
                  step <- min(c(h1,h2))/2
                  supp <- range(c(x,y))
                  p<- seq(from= supp[1], to =supp[2], by= step)        
                  dist<-sum(log((f(p)+me)/(g(p)+me))*f(p))*step
                  if(symmetrize)
                  {
                    dist <- (dist +  sum(log((g(p)+me)/(f(p)+me))*g(p))*step)/2
          
                  }
      		  return(dist)   


	      }
              else
              {
                  levs <- sort(unique(c(x, y)))
                  tabx <- table(factor(x, levels = levs)) / length(x)
                  taby <- table(factor(y, levels = levs)) / length(y)
                  dist<-sum(ifelse(tabx > 0, log(tabx / (taby + me)) * tabx, 0), na.rm = TRUE)/nbin
		  if(symmetrize)
                  {
		    dist <- (dist + sum(ifelse(taby > 0, log(taby / (tabx + me)) * taby, 0), na.rm = TRUE)/nbin)/2
		  }
              }
              return(dist)
          }
          ans<-rep(NA, n*(n-1)/2)
          ct <- 1
          for(i in 1:(n-1))
	      for(j in (i+1):n) {
		  if(!is.na(x[[i]]) && !is.na(x[[j]]))
                  	ans[ct] <- distfun(x[[i]], x[[j]])
		  else
		 	ans[ct]=NA
		  ct <- ct+1
	      }
          attributes(ans) <- list(Size = n, Labels = names(x),
                            Diag = diag, Upper = upper, methods =
                            "KLdist", class = "dist")
          
          ans
      })


# # # ## tentative "list" method for unequal sized samples (added by
# # # ## Deepayan Sarkar)
# # # 
# # # 
# # # setMethod("KLdist.matrix",
# # #           signature=signature("list"), 
# # #           function(x, 
# # #                    discretize = TRUE, nbin = 10,
# # #                    symmetrize = FALSE,
# # #                    diag = FALSE, upper=FALSE)
# # #       {
# # #           n <- length(x)
# # #           clist <- vector("list", length=n)
# # #           me <- .Machine$double.eps
# # #           
# # #           ##note: we combine x and y before binning, to make sure we span
# # #           ##   the range of the data, and we add machine epsilon to 
# # #           ##   protect against +/- Inf; this could use some work.
# # #           distfun <- function(x, y)
# # #           {
# # #               ## not clear what should be done if exactly one of x and y is a factor
# # #               if (discretize && !is.factor(x))
# # #               {
# # #                   breaks.x <- hist(c(x,y), breaks = nbin, plot = FALSE)$breaks
# # #                   temp1 <- table(cut(y, breaks.x, include.lowest = TRUE)) / length(y)
# # #                   ## temp1 <- temp1 + me
# # #                   temp2 <- table(cut(x, breaks.x, include.lowest = TRUE)) / length(x)
# # #                   ## temp2 <- temp2 + me
# # #                   sum( ifelse(temp2 > 0, log(temp2 / (temp1 + me)) * temp2, 0) , na.rm = TRUE)
# # #               }
# # #               else
# # #               {
# # #                   levs <- sort(unique(c(x, y)))
# # #                   tabx <- table(factor(x, levels = levs)) / length(x)
# # #                   taby <- table(factor(y, levels = levs)) / length(y)
# # #                   sum(ifelse(tabx > 0, log(tabx / (taby + me)) * tabx, 0), na.rm = TRUE)
# # #               }
# # #           }
# # #           ans <- matrix(NA, n, n)
# # #           for(i in seq_len(n))
# # #               for(j in seq_len(n))
# # #               {
# # #                   ans[i, j] <- distfun(x[[i]], x[[j]])
# # #               }
# # #           if(symmetrize) ans <- t(ans) + ans
# # #           if (!is.null(names(x)))
# # #               rownames(ans) <- colnames(ans) <- names(x)
# # #           ans
# # #       })
# # # 
