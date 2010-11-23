\name{sim.env.data}

\alias{sim.env.data}

\title{Simulates data for an environmental exposure}

\description{Generates data for a binary, quantitative-normal, or quantitative-uniform environmental determinant}

\usage{sim.env.data(num.obs = 20000, env.expo = 0, env.mean.lowlm = 3.3, 
env.stdev.uplm = 5, env.prev = 0.1, skewness = 0)}

\arguments{
  \item{num.obs}{Number of observations to simulate}
  \item{env.expo}{Model of the exposure: binary=0, quantitative-normal=1, quantitative-uniform=2}
  \item{env.mean.lowlm}{Mean under quantitative-normal model and lower limit under quantitative-uniform model}
  \item{env.stdev.uplm}{Standard deviation under quantitative-normal model and upper limit under quantitative-uniform model}
  \item{env.prev}{Prevalence of the environmental exposure}
  \item{skewness}{Determines skewness under quantitative-normal model; right-skewed distribution for positive values, left-skewed for negative value. The default is 0 (non skewed)}
}

\value{ A vector of continuous or binary values}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{make.obs.env}}}

\examples{
# Generate data for a binary enviromental exposure 
env.data <- sim.env.data(1000, 0, 3.3, 5, 0.1, 0)

}

\keyword{ utilities }
