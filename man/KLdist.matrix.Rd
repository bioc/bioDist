\name{KLdist.matriX}
\alias{KLdist.matrix}
\alias{KLdist.matrix,matrix-method}
\alias{KLdist.matrix,list-method}
\alias{KLdist.matrix,eSet-method}
\title{Discrete version of Kullback-Leibler Distance (KLD)}
\description{
Calculate the KLD by binning continuous data.

KL distance is calculated using the formula
     
\deqn{KLD(f_1(x),f_2(x)) = \sum_{i=1}^N{ 
       f_1(x_i)*\log\frac{f_1(x_i)}{f_2(x_i)}}}
}
                            

\usage{
KLdist.matrix(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or a list or an object of a class that extends eSet.
 If x is an object of a class derived from eSet (ExpressionSet,SnpSet etc),
then the values returned by the \code{exprs} function are used.}
  \item{\dots}{arguments passed to \code{KLdist.matrix}:
    \describe{
      \item{\code{gridsize}}{the number of grid points used to select the 
             optimal bin width of the histogram used to estimate density. 
             If no value is supplied, the grid size is calculated internally; 
             default is \code{NULL}.} 
      \item{\code{symmetrize}}{if \code{TRUE}, then symmetrize; the
         default is \code{FALSE}.}
      \item{\code{diag}}{if \code{TRUE}, then the diagonal of the distance 
           matrix will be displayed; the default is \code{FALSE}.} 
      \item{\code{upper}}{if \code{TRUE}, then the upper triangle of the 
         distance matrix will be displayed; default is \code{FALSE}.} 
      \item{\code{sample}}{for eSet methods: if \code{TRUE}, then the 
        distances are computed between samples, otherwise, between features; 
          the default is \code{TRUE}.}
    }
  }
}
\details{
  The data are binned, and then the KL distance between the two discrete
  distributions is computed and used. The distance is computed between
  rows of the input matrix (except if the input is an object of a class
  that extends eSet  and \code{sample} is \code{TRUE}.

  The presumption is that all samples have the same number of observations.
  The list method is meant for use when samples sizes are unequal.
}
\value{
An object of class \code{dist} is returned.
}
\author{Beiying Ding}
\seealso{\code{\link{cor.dist}}, \code{\link{spearman.dist}}, 
  \code{\link{tau.dist}},\code{\link{euc}},
  \code{\link{man}},\code{\link{KLD.matrix}},\code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(100), nrow = 5)
 KLdist.matrix(x, symmetrize = TRUE) 
}
\keyword{manip}
