\name{KLdist.matrix}
\alias{KLdist.matrix}
\alias{KLdist.matrix,matrix-method}
\alias{KLdist.matrix,list-method}
\alias{KLdist.matrix,ExpressionSet-method}
\title{Discrete version of Kullback-Leibler Distance (KLD)}
\description{
Calculate the KLD by binning continuous data.

KL distance is calculated using the formula
     
\deqn{KLD(f_1(x),f_2(x)) = \frac{1}{N}\sum_{i=1}^N{ 
       \Delta_i*f_1(x_i)*\log\frac{f_1(x_i)}{f_2(x_i)}}}
}
                            

\usage{
KLdist.matrix(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or ExpressionSet; if x is an ExpressionSet,then
    the values returned by the \code{exprs} function are used.}
  \item{\dots}{arguments passed to \code{KLdist.matrix}:
    \itemize{
      \item{nbin}{number of bins to calculate discrete probabilities; default
      is 10.} 
      \item{symmetrize}{if TRUE, then symmetrize; default is FALSE.}
      \item{diag}{if TRUE, then the diagonal of the distance matrix will be
      displayed; default is FALSE.} 
      \item{upper}{if TRUE, then the upper triangle of the distance matrix
      will be displayed; default is FALSE.} 
      \item{sample}{for ExpressionSet methods: if TRUE, then the distances are
      computed between samples, 
    otherwise, between genes.}
    }
  }
}
\details{
  The data are binned, and then the KL distance between the two discrete
  distributions is computed and used. The distance is computed between
  rows of the input matrix (except if the input is an ExpressionSet
  and \code{sample} is \cote{TRUE}.

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
 KLdist.matrix(x, symmetrize = TRUE, nbin = 3) 
}
\keyword{manip}
