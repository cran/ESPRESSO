\name{is.posdef}

\alias{is.posdef}

\title{Tells if a matrix is positive definite}

\description{
  Checks if any of the eigenvalues of the matrix is smaller than the set tolerance value.
}

\usage{is.posdef(matrix, tolerance = 1e-06)}

\arguments{
  \item{matrix}{Input matrix}
  \item{tolerance}{Constant}
}

\value{TRUE or FALSE}

\author{Amadou Gaye}

\seealso{\code{\link{make.posdef}}}

\examples{

# Example 1

# a positive definite matrix 
mat1 <- matrix(c(0.9999934,0.9999914,0.9999914,0.9999934),2,2) 

# check if the matrix is positive definite
is.posdef(mat1, 0.000001)

# Example 2

# a non positive definite matrix 
mat2 <- matrix(c(0.9999924,0.9999924,0.9999924,0.9999924),2,2) 

# check if the matrix is positive definite
is.posdef(mat2, 0.000001)

}

\keyword{utilities}

