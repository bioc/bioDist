\name{euc}
\alias{euc}
\title{ Euclidean distance} 
\description{
Calculate pairwise Euclidean distances for all columns of a matrix and save as a 'dist' object
}
\usage{
euc(x, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
   \item{diag}{if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{ if TRUE, upper triangle of the distance matrix will be displayed  } 
}
\details{
}
\value{
 Pairwise Euclidean distance object
}
\author{ Beiying Ding }

\seealso{ \code{\link{spearman.dist}}, \code{\link{tau.dist}},
 \code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},
 \code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(200), nrow=5)
 euc(x)
}
\keyword{  }
