\name{man}
\alias{man}
\title{ Manhattan distance} 
\description{
Calculate pairwise Manhattan distances for all columns of a matrix and save as a 'dist' object
}
\usage{
man(x, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
    \item{diag}{if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{ if TRUE, upper triangle of the distance matrix will be displayed  }
}
\details{
}
\value{
 Pairwise Manhattan distance object
}
\references{  }
\author{ Beiying Ding }
\note{ Requires the 'stats' package }


\seealso{ \code{\link{cor.dist}}, \code{link{spearman.dist}}, \code{\link{tau.dist}},\code{\link{euc}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}} }
\examples{
 x <- matrix(rnorm(200),ncol=5)
 man(x)
}
\keyword{  }
