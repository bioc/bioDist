\name{KLD.matrix}
\alias{KLD.matrix}
\title{ Continuous version of Kullback-Leibler Distance (KLD) }
\description{
Calculate KLD by estimating by smoothing log( f(x)/g(x))*f(x) and then integrate
}
\usage{
KLD.matrix(x, method = c("locfit", "density"), supp = c(-3, 3), subdivisions = 1000, diag = FALSE, upper = FALSE)
}
\arguments{
\item{x}{ n by p matrix }
  \item{method}{ use \code{locfit} or \code{density} to estimate integrand}
  \item{supp}{upper and lower limits of the integral }
  \item{subdivisions}{ subdivisions for the integration}
  \item{diag}{ if TRUE, diagonal of the distance matrix will be displayed }
  \item{upper}{if TRUE, upper triangle of the distance matrix will be displayed }
}
\details{
}
\value{
Pairwise KLD object
}
\references{  }
\author{ Beiying Ding, Vincent Carey }
\note{  }

\examples{
 x <- matrix(rnorm(100),ncol=5)
 KLD.matrix(x, method="locfit",supp=range(x))
}

\seealso{ \code{\link{cor.dist}}, \code{link{spearman.dist}}, \code{\link{tau.dist}},\code{\link{euc}},
\code{\link{man}}, \code{\link{KLdist.matrix}},\code{\link{mutualInfo}} } 
\keyword{}
