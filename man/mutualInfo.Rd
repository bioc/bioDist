\name{mutualInfo}
\alias{mutualInfo}
\alias{MIdist}
\title{Mutual Information} 
\description{
Calculate mutual information via binning
}
\usage{
mutualInfo(x, nbin = 10,  diag = FALSE, upper = FALSE)
MIdist(x, nbin = 10,  diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ n by p matrix }
  \item{nbin}{number of bins to calculate discrete probabilities }
  \item{diag}{ if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{if TRUE, upper triangle of the distance matrix will be displayed }
}
\details{
 For \code{mutualInfo} each row of \code{x} the data are divided into 
 \code{nbins} groups and then the mutual information is computed, treating 
 the data as if they were discrete.

 For \code{MIdist} we use the transformation proposed by Joe (1989),
\eqn{\delta^* = (1 - \exp(-2 \delta))^{1/2}}{delta* = (1 - exp(-2 delta))^.5}
where \eqn{\delta}{delta} is the mutual information. The \code{MIdist} is 
then \eqn{1 = \delta^*}{1-delta*}. Joe argues that this measure is then 
similar to Kendall's tau, \code{\link{tau.dist}}.
}
\value{
 An object of class \code{dist} which contains the pairwise distances.
}
\ref{H. Joe, Relative Entropy Measures of Multivariate Dependence,
JASA, 1989, 157-164.}
\author{ Robert Gentleman }

\seealso{ \code{\link{dist}}, \code{\link{KLdist.matrix}},
  \code{\link{cor.dist}}, \code{\link{KLD.matrix}} }
\examples{
 x <- matrix(rnorm(100),nrow=5)
 mutualInfo(x, nbin=3)
}
\keyword{manip}
