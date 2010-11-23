\name{skew.rnorm}

\alias{skew.rnorm}

\title{Allows to generate right or left-skewed normal distribution}

\description{This function allows one to set a level of skewness for normally distributed data.}

\usage{skew.rnorm(num.obs = 20000, mean = 0, sd = 1, skewness = 0)}

\arguments{
  \item{num.obs}{Number of observations to simulate}
  \item{mean}{Mean}
  \item{sd}{Standard deviation}
  \item{skewness}{Determines the direction and level of skewness; right-skewed for positive value, left-skewed for negative value. The default is 0 (non skewed)}
}

\value{A numerical vector}

\references{Azzalini, A. 1985, \code{A class of distributions which includes the normal ones.}, Scandinavian Journal of Statistics, \bold{vol. 12}, pp171-178.}

\author{Amadou Gaye}

\seealso{\code{\link{sim.env.data}}}

\examples{some.data <- skew.rnorm(20000, 0, 1, 0)}

\keyword{utilities}

