\name{make.posdef}

\alias{make.posdef}

\title{Turns a matrix into a positive definite one}

\description{
   Computes the nearest positive definite matrix of a real symmetric matrix.
}

\usage{make.posdef(matrix, tolerance = 1e-06)}

\arguments{
  \item{matrix}{Input matrix}
  \item{tolerance}{Constant}
}

\value{A positive-definite matrix}

\references{
N.J. Higham, 1988 \code{Computing a nearest symmetric positive semidefinite matrix}, Linear Algebra Appl. \bold{vol. 103}, pp.103 118
}

\author{Amadou Gaye}

\seealso{\code{\link{is.posdef}}}

\examples{

# a non positive definite matrix 
mat <- matrix(c(0.9999924,0.9999924,0.9999924,0.9999924),2,2)

# make the matrix positive definite  
mat.new <- make.posdef(mat, 0.000001)

}

\keyword{utilities}

