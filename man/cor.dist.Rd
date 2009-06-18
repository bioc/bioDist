\name{cor.dist}
\alias{cor.dist}
\alias{cor.dist,matrix-method}
\alias{cor.dist,eSet-method}
\title{Pearson correlational distance}
\description{
 Calculate pairwise Pearson correlational distances, i.e. 1-COR or 1-|COR|,
 and saves as a 'dist' object
}
\usage{
cor.dist(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or ExpressionSet; if x is an ExpressionSet, then the 
  function uses its 'exprs' slot.}
  \item{\dots}{arguments passed to \code{cor.dist}:
    \itemize{
      \item{abs}{if TRUE, then 1-|COR| else 1-COR, default is TRUE.}
      \item{diag}{if TRUE, then the diagonal of the distance matrix will be
      displayed, default is FALSE.} 
      \item{upper}{if TRUE, then the upper triangle of the distance matrix
      will be displayed, default is FALSE.} 
      \item{sample}{for objects of classes that extend eSet: if TRUE, then
	  distances are computed between samples(columns) , otherwise, they are
      computed between features(rows).}
    }
  }
}
\details{
The \code{cor} function is used to compute the pairwise distances between rows 
of an input matrix, except if the input is an object of a class that extends 
eSet and \code{sample} is \code{TRUE}.
}
\value{
  Pairwise Pearson correlational distance object
}
\author{Beiying Ding} 
\seealso{\code{\link{spearman.dist}}, \code{\link{tau.dist}},\code{\link{euc}},
\code{\link{man}},  \code{\link{KLdist.matrix}}, \code{\link{KLD.matrix}},
\code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(200), nrow = 5)
 cor.dist(x)
}
