\name{gen.params}

\alias{gen.params}

\docType{data}

\title{Parameters to simulate genetic data}

\description{
    A table of scenarios (rows) and parameters (columns).
}

\usage{data(gen.params)}

\format{
  A data frame with 24 observations on the following 14 variables.
  \describe{
    \item{\code{scenario.id}}{Scenario number; each row stores parameters for one scenario}
    \item{\code{model.geno1}}{Genetic model of the first variant: 0 for binary and 1 for additive}
    \item{\code{model.geno2}}{Genetic model of the second variant: 0 for binary and 1 for additive}
    \item{\code{MAF.geno1}}{Minor allele frequency of genetic variant 1}
    \item{\code{MAF.geno2}}{Minor allele frequency of genetic variant 2}
    \item{\code{or.geno1}}{Odds-ratio of genetic variant 1}
    \item{\code{geno1.efkt}}{effect of genetic variant 2}
    \item{\code{or.geno2}}{Odds-ratio of genetic variant 1}
    \item{\code{geno2.efkt}}{effect of genetic variant 2}
    \item{LD}{Sets independence or LD between the two genetic variants: 0 for independence and 1 for LD}
    \item{R.target}{Correlation coefficient required if the alleles of the two genetic variants are in LD}
    \item{\code{display}}{If TRUE, a summary is printed on screen}
    \item{\code{sensitivity.geno}}{Sensitivity of the assessment of genotype}
    \item{\code{specificity.geno}}{Specificity of the assessment of genotype}
  }
}

\examples{data(gen.params)}

\keyword{datasets}
