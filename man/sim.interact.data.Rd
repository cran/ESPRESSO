\name{sim.interact.data}

\alias{sim.interact.data}

\title{Generates data for the interaction term}

\description{
  Computes the interaction term for the pre-specified interaction model. 
}

\usage{sim.interact.data(geno1.U, geno2.U, env1.U, env2.U, is.interaction = 1)}

\arguments{
  \item{geno1.U}{Genotype data for genetic variant 1}
  \item{geno2.U}{Genotype data for genetic variant 2}
  \item{env1.U}{Exposure data for environment 1}
  \item{env2.U}{Exposure data for environment 2}
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
}

\value{A numerical vector}

\author{Amadou Gaye, Paul Burton}

\examples{

# genotypic data
geno1 <-sim.geno.data(10000, 0.1, 0)
geno2 <-sim.geno.data(10000, 0.2, 0)

# Environmental exposure data
env1 <- sim.env.data(1000, 0, 3.3, 5, 0.1, 0)
env2 <- sim.env.data(1000, 0, 3.3, 5, 0.2, 0)

# interaction data for a gene-environment interaction model
int.data <- sim.interact.data(geno1, geno2, env1, env2, 1)
}

\keyword{utilities}

