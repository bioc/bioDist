\name{euc}
\alias{euc}
\alias{euc,matrix-method}
\alias{euc,ExpressionSet-method}
\title{Euclidean distance} 
\description{
Calculate pairwise Euclidean distances for all columns of a matrix and save as a 'dist' object
}
\usage{
euc(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or ExpressionSet; if x is an ExpressionSet, then the
  function uses its 'exprs' slot.} 
  \item{\dots}{arguments passed to \code{euc}:
    \itemize{
       \item{\code{diag}}{if TRUE, then the diagonal of the distance matrix will
       be displayed; default is FALSE.} 
       \item{\code{upper}}{if TRUE, then the upper triangle of the distance
       matrix will be displayed; default is FALSE.}
    }
  }
}
\details{
}
\value{
 Pairwise Euclidean distance object
}
\author{Beiying Ding}
\seealso{\code{\link{spearman.dist}}, \code{\link{tau.dist}},
  \code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},
  \code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(200), nrow = 5)
 euc(x)
}
