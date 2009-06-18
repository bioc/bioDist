\name{man}
\alias{man}
\alias{man,matrix-method}
\alias{man,eSet-method}
\title{Manhattan distance} 
\description{
Calculate pairwise Manhattan distances and saves as a \code{dist} object.
}
\usage{
man(x, \dots)
}
\arguments{
  	\item{x}{n by p matrix or an object of class that extends eSet. If x is an object of class that extends eSet, (eg  ExpressionSet) then the function uses its 'exprs' slot.}
  	\item{\dots}{arguments passed to \code{man}:
    \itemize{
    \item{diag}{if \code{TRUE}, then the diagonal of the distance matrix will 
       be displayed; default is FALSE.}
    \item{upper}{if \code{TRUE}, then the upper triangle of the distance
      matrix will be displayed; default is FALSE.}
    }
  }
}
\details{
 This is just an interface to \code{dist} with the right parameters set.
}
\value{
 An instance of the \code{dist} class with the pairwise Manhattan distances
  between the rows of \code{x} in case of a matrix or between the features 
(rows) in case of a class that extends eSet.
}
\author{Beiying Ding}
\seealso{\code{\link{cor.dist}}, \code{\link{spearman.dist}}, 
  \code{\link{tau.dist}},\code{\link{euc}}, \code{\link{KLdist.matrix}},
  \code{\link{KLD.matrix}},\code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(200), nrow = 5)
 man(x)
}
\keyword{manip}
