\name{mutualInfo}
\alias{mutualInfo}
\title{Mutual Information} 
\description{
Calculate mutual information via binning
}
\usage{
mutualInfo(x, nbin = 10,  diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
  \item{nbin}{number of bins to calculate discrete probabilities }
  \item{diag}{ if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{if TRUE, upper triangle of the distance matrix will be displayed }
}
\details{
}
\value{
Pairwise mutual information distance object
}
\references{  }
\author{ Robert Gentleman }
\note{  }


\seealso{ \code{\link{cor.dist}}, \code{link{spearman.dist}}, \code{\link{tau.dist}},\code{\link{euc}},
\code{\link{man}}, \code{\link{KLdist.matrix}},\code{\link{KLD.matrix}} }    
\examples{
 x <- matrix(rnorm(100),ncol=5)
 mutualInfo(x, nbin=3)
}
\keyword{}
