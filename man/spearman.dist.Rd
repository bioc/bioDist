\name{spearman.dist}
\alias{spearman.dist}
\title{ Spearman correlational distance }
\description{
 Calculate pairwise Spearman correlational distances, i.e. 1-SPEAR or 1-|SPEAR|, for all columns of a matrix and save as a 'dist' object
}
\usage{
spearman.dist(x, abs = TRUE, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
  \item{abs}{ if TRUE, then 1-|SPEAR| else 1-SPEAR }
  \item{diag}{if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{ if TRUE, upper triangle of the distance matrix will be displayed  }
}
\details{
}
\value{
  Pairwise Spearman correlational distance object
}
\references{  }
\author{Beiying Ding} 
\note{  }

\seealso{ 
\code{\link{cor.dist}}, \code{\link{tau.dist}},\code{\link{euc}},\code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}}} 
\examples{
 x <- matrix(rnorm(200),ncol=5)
 spearman.dist(x)
}

\keyword{ }
