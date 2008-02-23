\name{man}
\alias{man}
\alias{man,matrix-method}
\alias{man,ExpressionSet-method}
\title{Manhattan distance} 
\description{
Calculate pairwise Manhattan distances for all rows of a matrix and 
save as a \code{dist} object.
}
\usage{
man(x, \dots)
}
\arguments{
  \item{x}{a matrix or an ExpressionSet; if x is an ExpressionSet, then the function uses its 'exprs' slot.}
  \item{\dots}{arguments passed to \code{man}:
    \item{diag}{if \code{TRUE}, then the diagonal of the distance matrix will 
     be displayed; default is FALSE.}
    \item{upper}{if \code{TRUE}, then the upper triangle of the distance matrix will 
     be displayed; default is FALSE.}
}
}
\details{
 This is just an interface to \code{dist} with the right parameters set.
}
\value{
 An instance of the \code{dist} class with the pairwise Manhattan distances
  between the rows of \code{x}.
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
