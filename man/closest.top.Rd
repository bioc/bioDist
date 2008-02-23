\name{closest.top}
\alias{closest.top}
\title{Find the closest genes.}
\description{
Find the closest genes to the supplied target gene based on the
supplied distances.
}
\usage{
closest.top(x, dist.mat, top)
}
\arguments{
  \item{x}{the name of the gene (feature) to use.}
  \item{dist.mat}{either a dist object or a matrix of distances.}
  \item{top}{the number of closest genes desired.}
}
\details{
The feature named \code{x} must be in the supplied distances. If so, then 
the \code{top} closest other features are returned.
}
\value{
 A vector of names of the \code{top} closest features.
}
\author{Beiying Ding}
\seealso{\code{\link{cor.dist}}, \code{\link{spearman.dist}}, \code{\link{tau.dist}},\code{\link{euc}},
  \code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}}}
\examples{
 data(sample.ExpressionSet)
 sE <- sample.ExpressionSet[1:100,]
 d1 <- KLdist.matrix(sE, sample = FALSE)
 closest.top(featureNames(sE)[1], d1, 5)
}
