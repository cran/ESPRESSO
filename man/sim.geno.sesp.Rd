\name{sim.geno.sesp}

\alias{sim.geno.sesp}

\title{Simulates sensitivity and specificity for genotype assessment}

\description{
  Generates the sensitivity and specificity values to assess the genotype of an unobserved variant based on the knowledge of an observed variant in LD with the unobserved one. 
} 

\usage{sim.geno.sesp(seed.val = 333333, prevalence.exp = 0.5, R2.target = 0.8)}

\arguments{
  \item{seed.val}{Seed value}
  \item{prevalence.exp}{prevalence of the observed genetic variant}
  \item{R2.target}{Measure of LD between the observed and the unobserved variant}
}

\value{
  A vector of two values:
  \item{sensitivity.mid}{Simulated sensitivity}
  \item{specificity.mid}{Simulated specificity}
}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{sim.env.sesp}}}

\examples{
   # simulate sensitivity and specificity
   geno.sens.spec <- sim.geno.sesp(333333,0.5,0.8)
}

\keyword{utilities}

