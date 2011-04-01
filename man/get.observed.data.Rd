\name{get.observed.data}

\alias{get.observed.data}

\title{Generates exposure data with with some error}

\description{Uses functions \code{make.obs.geno} and \code{make.obs.env} to generate effect data with a set level of error.}

\usage{
get.observed.data(is.interaction = 0, true.data, geno.error = c(0.05, 0.05), 
is.add = c(0, 0), MAF = c(0.1, 0.1), env.expo = c(0, 0), env.prev = c(0.1, 0.1), 
env.error = c(0.15, 0.15), reliability.env=c(0.8,0.8))
}

\arguments{
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{true.data}{Input table of simulated data considered as \code{true} data}
  \item{geno.error}{Misclassification rates in the assessment of genotypes}
  \item{is.add}{Genetic models of the two variants: 0 for binary model and 1 for additive model}
  \item{MAF}{Minor Allele frequencies}
  \item{env.expo}{Model of the exposure: binary=1, quantitative-normal=1 or quantitative-uniform=2}
  \item{env.prev}{Prevalence of the two environmental exposures}
  \item{env.error}{Misclassification rates in environmental exposures assessment: \code{1-sensitivity} and \code{1-specificity}}
  \item{reliability.env}{Reliability of the assessment of quantitative exposures}
}

\value{A matrix containing 11 variables}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{make.obs.geno}}, \code{\link{make.obs.env}}}

\examples{
# load the 'true' data 
data(init.data)

# adds some error to the 'true' exposure data and generate 'observed' data
obs.data <- get.observed.data(0, init.data, c(0.05, 0.05), c(0, 0), c(0.1, 0.1), 
c(0, 0), c(0.1, 0.1), c(0.15, 0.15), c(0.8, 0.8))

}

\keyword{utilities}

