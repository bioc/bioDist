\name{cor.dist}
\alias{cor.dist}
\alias{cor.dist,matrix-method}
\alias{cor.dist,exprSet-method}
\alias{cor.dist,ExpressionSet-method}
\title{ Pearson correlational distance }
\description{
 Calculate pairwise Pearson correlational distances, i.e. 1-COR or 1-|COR|, for all columns of a matrix and save as a 'dist' object
}
\usage{
cor.dist(x, \dots)
}
\arguments{
  \item{x}{ n by p matrix or exprSet; if x is an exprSet, function workes against its 'exprs' slot}
  \item{\dots}{arguments passed to \code{cor.dist}:
  \item{abs}{ if TRUE, then 1-|COR| else 1-COR, default is TRUE }
  \item{diag}{if TRUE, diagonal of the distance matrix will be displayed, default is FALSE }
  \item{upper}{ if TRUE, upper triangle of the distance matrix will be displayed, default is FALSE  }
}
}
\details{
}
\value{
  Pairwise Pearson correlational distance object
}
\author{Beiying Ding} 

\seealso{ \code{\link{spearman.dist}}, \code{\link{tau.dist}},\code{\link{euc}},\code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}} }
\examples{
 x <- matrix(rnorm(200), nrow=5)
 cor.dist(x)
}

\keyword{ }
