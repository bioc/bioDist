\name{KLdist.matrix}
\alias{KLdist.matrix}
\title{ Discrete version of Kullback-Leibler Distance (KLD)}
\description{
Calculate KLD by binning
}
\usage{
KLdist.matrix(x, nbin = 10, symmetrize = TRUE, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
  \item{nbin}{number of bins to calculate discrete probabilities } 
\item{symmetrize}{if TRUE, then symmetrize }
  \item{diag}{ if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{if TRUE, upper triangle of the distance matrix will be displayed }
}
\details{
}

\value{
Pairwsie KLD object
}
\author{ Beiying Ding}

\seealso{ \code{\link{cor.dist}}, \code{\link{spearman.dist}}, 
\code{\link{tau.dist}},\code{\link{euc}},
\code{\link{man}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}} } 
\examples{
 x <- matrix(rnorm(100), nrow=5)
 KLdist.matrix(x,symmetrize=TRUE, nbin=3) 
}
\keyword{manip}
