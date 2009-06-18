\name{KLD.matrix}
\alias{KLD.matrix}
\alias{KLD.matrix,matrix-method}
\alias{KLD.matrix,list-method}
\alias{KLD.matrix,eSet-method}
\title{Continuous version of Kullback-Leibler Distance (KLD)}
\description{
Calculate KLD by estimating by smoothing 
 \eqn{\log(f(x)/g(x))*f(x)}{\log(f(x)/g(x))*f(x)} and then integrating.
}
\usage{
KLD.matrix(x, \dots) 
}
\arguments{
  \item{x}{n by p matrix or list or an object of a class that extends eSet; if x is an an object of 
  a class that extends eSet (eg ExpressionSet), then the function works against its 'exprs' slot.}
  \item{\dots}{arguments passed to \code{KLD.matrix}:
    \itemize{
      \item{method}{use \code{locfit} or \code{density} to estimate integrand;
      default is c("locfit", "density")(i.e. both methods).} 
      \item{supp}{upper and lower limits of the integral; default is NULL in
      which case the limits of the integral are calculated from the range of 
      the data.}
      \item{subdivisions}{subdivisions for the integration; default is 1000.}
      \item{diag}{if TRUE, then the diagonal of the distance matrix will be
      displayed; default is FALSE.} 
      \item{upper}{if TRUE, then the upper triangle of the distance matrix
      will be displayed; default is FALSE.} 
      \item{sample}{for ExpressionSet methods: if TRUE, then distances are
      computed between samples, otherwise, they are computed between genes.}
    }
  }
}
\details{
  The distance is computed between rows of the input matrix 
  (except if the input is an object of a class  that extends eSet 
   and \code{sample} is \code{TRUE}.

  The presumption is that all samples have the same number of observations.
  The list method is meant for use when samples sizes are unequal.

}
\value{
An object of class \code{dist} with the pairwise, between rows, 
Kullback-Leibler distances.
}
\author{Beiying Ding, Vincent Carey}
\seealso{ \code{\link{cor.dist}}, \code{\link{spearman.dist}}, 
  \code{\link{tau.dist}}, \code{\link{dist}}, 
  \code{\link{KLdist.matrix}}, \code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(100), nrow = 5)
 KLD.matrix(x, method = "locfit", supp = range(x))
}
\keyword{manip}
