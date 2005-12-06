\name{KLdist.matrix}
\alias{KLdist.matrix}
\alias{KLdist.matrix,matrix-method}
\alias{KLdist.matrix,exprSet-method}
\title{ Discrete version of Kullback-Leibler Distance (KLD)}
\description{
Calculate KLD by binning
}
\usage{
KLdist.matrix(x, \dots)
}
\arguments{
  \item{x}{ n by p matrix or exprSet; if x is an exprSet, function workes against its 'exprs' slot}
  \item{\dots}{arguments passed to \code{KLdist.matrix}:
  \item{nbin}{number of bins to calculate discrete probabilities, default is 10 } 
\item{symmetrize}{if TRUE, then symmetrize, default is FALSE }
  \item{diag}{ if TRUE, diagonal of the distance matrix will be displayed, default is FALSE }
  \item{upper}{if TRUE, upper triangle of the distance matrix will be displayed, default is FALSE }
}
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
