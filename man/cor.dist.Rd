\name{cor.dist}
\alias{cor.dist}
\title{ Pearson correlational distance }
\description{
 Calculate pairwise Pearson correlational distances, i.e. 1-COR or 1-|COR|, for all columns of a matrix and save as a 'dist' object
}
\usage{
cor.dist(x, abs = TRUE, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
  \item{abs}{ if TRUE, then 1-|COR| else 1-COR }
  \item{diag}{if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{ if TRUE, upper triangle of the distance matrix will be displayed  }
}
\details{
}
\value{
  Pairwise Pearson correlational distance object
}
\references{  }
\author{Beiying Ding} 
\note{  }

\seealso{ \code{\link{spearman.dist}}, \code{\link{tau.dist}},\code{\link{euc}},\code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}} }
\examples{
 x <- matrix(rnorm(200),ncol=5)
 cor.dist(x)
}

\keyword{ }
