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
 For each row of \code{x} the data are divided into \code{nbins} groups and 
 then the mutual information is computed, treating the data as if they were
 discrete.
}
\value{
 An object of class \code{dist} which contains the pairwise distances.
}
\author{ Robert Gentleman }


\seealso{ \code{\link{dist}}, \code{\link{KLdist.matrix}},
  \code{\link{cor.dist}}, \code{\link{KLD.matrix}} }
\examples{
 x <- matrix(rnorm(100),nrow=5)
 mutualInfo(x, nbin=3)
}
\keyword{}
