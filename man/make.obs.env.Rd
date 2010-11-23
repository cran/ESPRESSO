\name{make.obs.env}

\alias{make.obs.env}

\title{Adds some error to environmental exposure data}

\description{
  Adds a set level of error to simulated binary or quantitative data (the true data) to obtain data with a larger variance (the observed data). The level of error is determined by the misclassification rates in binary data and by the set level of variance in the quantitative data.
}

\usage{
make.obs.env(env.data, env.expo = 0, env.prev = 0.1, env.error = c(0.15, 0.15), 
stdev.error = 0.18, env.obs.var = 0.04)
}

\arguments{
  \item{env.data}{ A vector of environmental measures that represents the true data }
  \item{env.expo}{ Model of the exposure: binary=1, quantitative-normal=1 or quantitative-uniform=2 }
  \item{env.prev}{ Prevalence of the environmental determinant }
  \item{env.error}{ Misclassification rates: \code{1-sensitivity} and \code{1-specificity} }
  \item{stdev.error}{ Standard deviation of the error in the observed data }
  \item{env.obs.var}{ Variance of the observed data }
}

\value{
  A dataframe of two columns:
  \item{environ.orig}{input data (true data)}
  \item{environ.new}{observed data}
}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{sim.env.data}}, \code{\link{misclassify}}}

\examples{

# load a dataframe that contains environmental exposure data
data(init.data)

# get observed data by adding some error to the initial data
true.env.data <- init.data$env1.U
x <- make.obs.env(true.env.data)
observed.env.data <- x$environ.new

}

\keyword{utilities}

