\name{sim.pheno.bin}

\alias{sim.pheno.bin}

\title{Simulates binary outcome data}

\description{
  Uses the effects data and the odds-ratios of genetic, environmental and eventually interaction determinants to contruct a linear predictor (LP). The probability of disease, obtained through logictic transformation of the LP, is used to generate a bionomially distributed outcome (the true phenotypes). The vector of true phenotypes is then randomly misclassified to obtained the observed phenotypes. The level of misclassification is given by the sensitivity and specificity of the phenotype assessment.
}

\usage{
sim.pheno.bin(num.obs = 20000, is.interaction = 0, disease.prev = 0.1, geno1.U, 
geno2.U, env1.U, env2.U, int.U, subject.effect, or.geno = c(1.5, 1.5), 
or.env = c(1.5, 1.5), or.int = 1.8, pheno.error = c(0, 0))
}


\arguments{
  \item{num.obs}{Number of observations to simulate}
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{disease.prev}{Prevalence of the binary outcome}
  \item{geno1.U}{Genotype data for genetic variant 1}
  \item{geno2.U}{Genotype data for genetic variant 2}
  \item{env1.U}{Exposure data for environment 1}
  \item{env2.U}{Exposure data for environment 2}
  \item{int.U}{Interaction effect data}
  \item{subject.effect}{Subject effect data, reflects the heterogeneity in baseline disease risk}
  \item{or.geno}{Odds ratios of the two genetic variants}
  \item{or.env}{Odds ratios of the two environments}
  \item{or.int}{Odds ratio of the interaction}
  \item{pheno.error}{Misclassification rates in phenotype assessment}
}

\value{
  A dataframe of two columns:
  \item{pheno.original}{true phenotypes}
  \item{environ.new}{observed phenotypes}
}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{misclassify}}}

\examples{

# simulate genotype data for 2 genetic variants
geno1.data <- sim.geno.data(10000, 0.1, 0)
geno1 <- geno1.data$geno.U
geno2.data <- sim.geno.data(10000, 0.2, 0)
geno2 <- geno2.data$geno.U

# simulate environmental measures for two environments 
env1 <- sim.env.data(10000, 0, 3.3, 5, 0.1, 0)
env2 <- sim.env.data(10000, 0, 3.8, 5, 0.1, 0)

# generate interaction effect data
int <- sim.interact.data(geno1, geno2, env1, env2, 1)

# simulate subject effect data
subject.effect <- sim.subject.data(10000, 12.36)

# generate phenotypes
pheno.data <- sim.pheno.bin(10000, 1, 0.1, geno1, geno2, env1, env2, int, 
subject.effect, c(1.5,1.5), c(1.5,1.5), 1.8, c(0,0))

true.pheno <- pheno.data$pheno.original
observed.pheno <- pheno.data$pheno.U

}

\keyword{utilities}

