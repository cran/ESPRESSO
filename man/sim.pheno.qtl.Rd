\name{sim.pheno.qtl}

\alias{sim.pheno.qtl}

\title{Simulates continuous outcome data}

\description{
   Uses the effects data of the determinants to contruct a linear predictor (LP). The outcome is normally distributed variable generated with a mean equal to LP and a standard deviation of 1. Some error is then added to the simulated outcome to obtained the observed outcome.
}

\usage{
sim.pheno.qtl(num.subjects = 10000, is.interaction = 0, geno1.U, geno2.U, 
env1.U, env2.U, int.U, geno.efkt = c(0.25, 0.25), env.efkt = c(0.25, 0.25), 
int.efkt = 0.5, reliability.pheno)
}

\arguments{
  \item{num.subjects}{Number of subjects to simulate}
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{geno1.U}{Genotypes for genetic variant 1}
  \item{geno2.U}{Genotypes for genetic variant 2}
  \item{env1.U}{Exposure data for environment 1}
  \item{env2.U}{Exposure data for environment 2}
  \item{int.U}{Interaction effect data}
  \item{geno.efkt}{Effects of the genetic variants}
  \item{env.efkt}{Effects of the environmental determinants}
  \item{int.efkt}{Effect of the interaction term}
  \item{reliability.pheno}{Reliability of the phenotype assessment}
}

\value{
  A dataframe of two columns:
  \item{pheno.original}{true phenotypes}
  \item{environ.new}{observed phenotypes}
}

\author{Amadou Gaye, Paul Burton}

\examples{

# simulate genotype data for 2 genes
geno1.data <- sim.geno.data(10000, 0.1, 0)
geno1 <- geno1.data$geno.U
geno2.data <- sim.geno.data(10000, 0.2, 0)
geno2 <- geno2.data$geno.U

# simulate environmental measures for two environments 
env1 <- sim.env.data(10000, 0, 3.3, 5, 0.1, 0)
env2 <- sim.env.data(10000, 0, 3.8, 5, 0.1, 0)

# generate interaction effect data
int <- sim.interact.data(geno1, geno2, env1, env2, 1)

# generate phenotypes
pheno.data <- sim.pheno.qtl(10000, 1, geno1, geno2, env1, env2, int, 
c(0.25, 0.25), c(0.25, 0.25), 0.5, 0.9)

true.pheno <- pheno.data$pheno.original
observed.pheno <- pheno.data$pheno.U

}

\keyword{utilities}

