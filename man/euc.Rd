\name{euc}
\alias{euc}
\alias{euc,matrix-method}
\alias{euc,eSet-method}
\title{Euclidean distance} 
\description{
Calculate pairwise Euclidean distances and saves the result  as a 'dist' object
}
\usage{
euc(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or an object of a class that extends eSet; if x is a
  matrix, pairwise distances are calculated between the rows of a matrix. 
  If x is an object of a class that extends eSet, the method makes use of the
  'exprs' method and  pairwise distances are calculated
  between samples(columns) if \code{sample} is \code{TRUE}}   
  \item{\dots}{arguments passed to \code{euc}:
    \itemize{
       \item{\code{diag}}{if TRUE, then the diagonal of the distance matrix will
       be displayed; default is FALSE.} 
       \item{\code{upper}}{if TRUE, then the upper triangle of the distance
       matrix will be displayed; default is FALSE.}
       \item{\code{sample}}{For objects of classes that extends eSet, pairwise
       distances are calculated between samples(columns) if \code{sample} is 
	   \code{TRUE} ; default value is TRUE}
       }

   }
}
\details{
The method calculates pairwise euclidean distances, assuming that all samples 
have the same number of observations 
}
\value{
 An object of class \code{dist} with the pairwise Euclidean distance between rows
 except in case of objects of class that extend eSet when \code{sample} is 
\code{TRUE}
}
\author{Beiying Ding}
\seealso{\code{\link{spearman.dist}}, \code{\link{tau.dist}},
  \code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},
  \code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(200), nrow = 5)
 euc(x)
}
