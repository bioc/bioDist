mutualInfo <- function(x, nbin=10, diag=FALSE, upper=FALSE)
{
   x <- as.matrix(x)
   nc <- ncol(x)
   nr <- nrow(x)
   clist <- vector("list", length=nr)
   for(i in 1:nr)
       clist[[i]] <- cut(x[i,], breaks=nbin)

   ppfun <- function(pp) {pp<-pp[pp>0]; -sum(pp*log(pp ))}
   appfun <- function(x,y) {ppfun(table(x)/nc)+ppfun(table(y)/nc) -
                                 ppfun(c(table(x, y)/nc))}

   rvec<-rep(NA, nr*(nr-1)/2)
   ct <- 1
   for(i in 1:(nr-1))
       for(j in (i+1):nr) {
           rvec[ct] <- appfun(clist[[i]], clist[[j]])
           ct <- ct+1
   }   
   attributes(rvec) <- list(Size = nr, Labels = row.names(x),
                            Diag = diag, Upper = upper, methods =
                            "mutualInfo", class = "dist")
   rvec
}

MIdist = function(x, nbin=10, diag=FALSE, upper=FALSE) 
  1 - (1 - exp(-2*mutualInfo(x, nbin, diag, upper)))^.5
