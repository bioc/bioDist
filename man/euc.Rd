\name{euc}
\alias{euc}
\alias{euc,matrix-method}
\alias{euc,exprSet-method}
\alias{euc,ExpressionSet-method}
\title{ Euclidean distance} 
\description{
Calculate pairwise Euclidean distances for all columns of a matrix and save as a 'dist' object
}
\usage{
euc(x, \dots)
}
\arguments{
  \item{x}{n by p matrix or exprSet; if x is an exprSet, function workes against its 'exprs' slot}
  \item{\dots}{arguments passed to \code{euc}:
   \item{\code{diag}}{if TRUE, diagonal of the distance matrix will be displayed, default is FALSE }
  \item{\code{upper}}{ if TRUE, upper triangle of the distance matrix will be displayed, default is FALSE  } 
	}	
}
\details{
}
\value{
 Pairwise Euclidean distance object
}
\author{ Beiying Ding }

\seealso{ \code{\link{spearman.dist}}, \code{\link{tau.dist}},
 \code{\link{man}},\code{\link{KLdist.matrix}},\code{\link{KLD.matrix}},
 \code{\link{mutualInfo}}}
\examples{
 x <- matrix(rnorm(200), nrow=5)
 euc(x)
}
\keyword{  }
