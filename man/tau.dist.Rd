\name{tau.dist}
\alias{tau.dist}
\title{ Kendall's tau correlational distance }
\description{
 Calculate pairwise Kendall's tau correlational distances, i.e. 1-TAU or 1-|TAU|, for all columns of a matrix and save as a 'dist' object
}
\usage{
tau.dist(x, abs = TRUE, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
  \item{abs}{ if TRUE, then 1-|TAU| else 1-TAU }
  \item{diag}{if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{ if TRUE, upper triangle of the distance matrix will be displayed  }
}
\details{
}
\value{
  Pairwise Kendall's tau correlational distance object
}
\references{  }
\author{Beiying Ding} 
\note{  }

\seealso{ \code{\link{cor.dist}}, \code{\link{spearman.dist}},\code{\link{euc}},\code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}} }
\examples{
 x <- matrix(rnorm(200),ncol=5)
 tau.dist(x)
}

\keyword{ }
