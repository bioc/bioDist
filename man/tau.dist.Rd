\name{tau.dist}
\alias{tau.dist}
\alias{tau.dist,matrix-method}
\alias{tau.dist,exprSet-method}
\title{ Kendall's tau correlational distance }
\description{
 Calculate pairwise Kendall's tau correlational distances, i.e. 1-TAU or 1-|TAU|, 
 for all rows of the input matrix and return an instance of the 
 \code{dist} class.}
}
\usage{
tau.dist(x)
}
\synopsis{
tau.dist(x, \dots)
}
\arguments{
  \item{x}{ n by p matrix or exprSet; if x is an exprSet, function workes against its 'exprs' slot}
  \item{\dots}{arguments passed to \code{tau.dist}:
  \item{abs}{ if TRUE, then 1-|TAU| else 1-TAU, default is TRUE }
  \item{diag}{if TRUE, diagonal of the distance matrix will be displayed, default is FALSE }
  \item{upper}{ if TRUE, upper triangle of the distance matrix will be displayed, default is FALSE }
}
}
\details{
  Row-wise correlations are computed by calling the \code{cor} function 
 with the appropriate arguments.
}
\value{
  One minus the row-wise Kendall's tau correlations are returned as an 
 instance of the \code{dist} class.
}
\author{Beiying Ding} 

\seealso{ \code{\link{cor.dist}}, \code{\link{spearman.dist}},
\code{\link{euc}}, \code{\link{man}}, \code{\link{KLdist.matrix}},
\code{\link{KLD.matrix}}, \code{\link{mutualInfo}} }
\examples{
 x <- matrix(rnorm(200), nrow=5)
 tau.dist(x)
}

\keyword{manip}
