\name{sim.env.sesp}

\alias{sim.env.sesp}

\title{Simulates sensitivity and specificity for environmental exposure assessment}

\description{
  Generates environmental exposure status and tabulates truly exposed versus observed exposed to determine the number of false positives and false negatives which are then used to compute sensitivity and specificity.   
}

\usage{sim.env.sesp(seed.val = 333333, prevalence.exp = 0.5, reliability = 0.8)}

\arguments{
  \item{seed.val}{Seed value}
  \item{prevalence.exp}{Prevalence of the exposure}
  \item{reliability}{Reliability of the environmental exposure assessment}
}

\details{
  The function uses reliability as estimate of a standardized error (\code{1-reliability}) and generates exposure data with and without error. The prevalence of the exposure represents the probability in a quantile function to determine the threshold of truly exposed in the observed data (data without error) and in the true data (data with error). 
}

\value{
  A vector of two values:
  \item{sensitivity}{Simulated sensitivity of environmental exposure assessment}
  \item{specificity}{Simulated specificity of environmental exposure assessment}
}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{sim.geno.sesp}}}

\examples{
   # simulate sensitivity and specificity
   env.sens.spec <- sim.env.sesp(333, 0.2, 0.8)
}

\keyword{ utilities }
