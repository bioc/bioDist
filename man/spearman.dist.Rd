\name{spearman.dist}
\alias{spearman.dist}
\alias{spearman.dist,matrix-method}
\alias{spearman.dist,ExpressionSet-method}
\title{Spearman correlational distance}
\description{
Calculate pairwise Spearman correlational distances, i.e. 1-SPEAR or 
1-|SPEAR|, for all rows of a matrix and return a \code{dist} object.
}
\usage{
spearman.dist(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or ExpressionSet; if x is an ExpressionSet, then the function uses its 'exprs' slot.}
  \item{\dots}{arguments passed to \code{spearman.dist}:
  \item{abs}{if \code{TRUE}, then 1-|SPEAR| else 1-SPEAR; default is TRUE.}
  \item{diag}{if \code{TRUE}, then the diagonal of the distance matrix will be displayed; default is FALSE.}
  \item{upper}{if \code{TRUE}, then the upper triangle of the distance matrix will be displayed; default is FALSE.}
  \item{sample}{for the ExpressionSet method: if TRUE (the default), then distances are computed between samples.}
 }
}
\details{
We call \code{cor} with the appropriate arguments to compute the row-wise correlations. 
}
\value{
One minus the Spearman correlation, between rows of \code{x}, are returned, as 
an instance of the \code{dist} class.
}
\author{Beiying Ding} 
\seealso{
  \code{\link{cor.dist}}, \code{\link{tau.dist}}, \code{\link{euc}},
  \code{\link{man}}, \code{\link{KLdist.matrix}}, \code{\link{KLD.matrix}},
  \code{\link{mutualInfo}}, \code{\link{dist}}} 
\examples{
 x <- matrix(rnorm(200), nrow = 5)
 spearman.dist(x)
}
\keyword{manip}
