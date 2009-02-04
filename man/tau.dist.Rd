\name{tau.dist}
\alias{tau.dist}
\alias{tau.dist,matrix-method}
\alias{tau.dist,ExpressionSet-method}
\title{Kendall's tau correlational distance}
\description{
 Calculate pairwise Kendall's tau correlational distances, 
 i.e. 1-TAU or 1-|TAU|, 
 for all rows of the input matrix and return an instance of the 
 \code{dist} class.
}
\usage{
tau.dist(x)
}
\synopsis{
tau.dist(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or ExpressionSet; if x is an ExpressionSet, then the function uses its 'exprs' slot.}
  \item{\dots}{arguments passed to \code{tau.dist}:
    \itemize{
      \item{abs}{if TRUE, then 1-|TAU| else 1-TAU; default is TRUE.}
      \item{diag}{if TRUE, then the diagonal of the distance matrix will be
      displayed; default is FALSE.} 
      \item{upper}{if TRUE, then the upper triangle of the distance matrix
      will be displayed; default is FALSE.} 
      \item{sample}{for the ExpressionSet method: if TRUE (the default), then
      distances are computed between samples.} 
    }
  }
}
\details{
 Row-wise correlations are computed by calling the \code{cor} function 
 with the appropriate arguments.
}
\value{
 One minus the row-wise Kendall's tau correlations are returned as an 
 instance of the \code{dist} class. Note that this can be extremely 
 slow for large data sets.
}
\author{Beiying Ding} 
\seealso{\code{\link{cor.dist}}, \code{\link{spearman.dist}},
  \code{\link{euc}}, \code{\link{man}}, \code{\link{KLdist.matrix}},
  \code{\link{KLD.matrix}}, \code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(200), nrow = 5)
 tau.dist(x)
}
\keyword{manip}
