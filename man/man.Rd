\name{man}
\alias{man}
\title{ Manhattan distance} 
\description{
Calculate pairwise Manhattan distances for all rows of a matrix and 
save as a \code{dist} object.
}
\usage{
man(x, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ A matrix. }
    \item{diag}{If \code{TRUE}, the diagonal of the distance matrix will 
     be displayed }
  \item{upper}{If \code{TRUE}, the upper triangle of the distance matrix will 
     be displayed  }
}
\details{
 This is just an interface to \code{dist} with the right parameters set.
}
\value{
 An instance of the \code{dist} class with the pairwise Manhattan distances
  between the rows of \code{x}.
}
\author{ Beiying Ding }

\seealso{ \code{\link{cor.dist}}, \code{\link{spearman.dist}}, 
 \code{\link{tau.dist}},\code{\link{euc}}, \code{\link{KLdist.matrix}},
  \code{\link{KLD.matrix}},\code{\link{mutualInfo}} }
\examples{
 x <- matrix(rnorm(200), nrow = 5)
 man(x)
}
\keyword{manip}
