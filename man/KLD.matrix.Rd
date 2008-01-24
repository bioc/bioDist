\name{KLD.matrix}
\alias{KLD.matrix}
\alias{KLD.matrix,matrix-method}
\alias{KLD.matrix,list-method}
\alias{KLD.matrix,exprSet-method}
\alias{KLD.matrix,ExpressionSet-method}
\title{Continuous version of Kullback-Leibler Distance (KLD) }
\description{
Calculate KLD by estimating by smoothing 
 \eqn{\log( f(x)/g(x))*f(x)}{log( f(x)/g(x))*f(x)} and then integrating.
}
\usage{
KLD.matrix(x, \dots) 
}
\arguments{
  \item{x}{ n by p matrix or ExpressionSet; if x is an ExpressionSet, the function works against its 'exprs' slot}
  \item{\dots}{arguments passed to \code{KLD.matrix}:
  \item{method}{ use \code{locfit} or \code{density} to estimate integrand, default is c("locfit", "density")(i.e. both methods)}
  \item{supp}{upper and lower limits of the integral, default is c(-3, 3) }
  \item{subdivisions}{ subdivisions for the integration, default is 1000}
  \item{diag}{ if TRUE, diagonal of the distance matrix will be displayed, default is FALSE }
  \item{upper}{if TRUE, upper triangle of the distance matrix will be displayed, default is FALSE }
  \item{sample}{for ExpressionSet methods, if TRUE then distances are computed
   between samples, otherwise they are computed between genes.}
}
}
\details{
  The pairwise distances between the rows of \code{x} are computed.

  The list method is meant for use when samples sizes are unequal.
}
\value{
An object of class \code{dist} with the pairwise, between rows, 
Kullback-Leibler distances.
}
\author{ Beiying Ding, Vincent Carey }

\seealso{ \code{\link{cor.dist}}, \code{\link{spearman.dist}}, 
  \code{\link{tau.dist}}, \code{\link{dist}}, 
  \code{\link{KLdist.matrix}},\code{\link{mutualInfo}} } 
\examples{
 x <- matrix(rnorm(100),nrow=5)
 KLD.matrix(x, method="locfit",supp=range(x))
}

\keyword{manip}
