\name{spearman.dist}
\alias{spearman.dist}
\title{ Spearman correlational distance }
\description{
 Calculate pairwise Spearman correlational distances, i.e. 1-SPEAR or 
1-|SPEAR|, for all rows of a matrix and return a \code{dist} object.
}
\usage{
spearman.dist(x, abs = TRUE, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
  \item{abs}{ if \code{TRUE}, then 1-|SPEAR| else 1-SPEAR }
  \item{diag}{if \code{TRUE}, diagonal of the distance matrix will be displayed }
  \item{upper}{ if \code{TRUE}, upper triangle of the distance matrix will 
    be displayed  }
}
\details{
  We call \code{cor} with the appropriate arguments to compute the row-wise
correlations. 
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
 x <- matrix(rnorm(200), nrow=5)
 spearman.dist(x)
}

\keyword{manip}
