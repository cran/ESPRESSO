\name{init.data}

\alias{init.data}

\docType{data}

\title{Simulated genotypes, phenotypes and enviromental exposure data}

\description{
  A table of simulated \code{true} data. The number of rows represents the number of cases and controls (under binary outcome) or the number of subjects (under continuous outcome).
}

\usage{data(init.data)}

\format{
  A data frame with 10000 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{Scenario id}
    \item{\code{cc.U}}{Phenotypes}
    \item{\code{geno1.U}}{Genotype for the first genetic variant}
    \item{\code{geno2.U}}{Genotype for the second genetic variant}
    \item{\code{allele.A1}}{Alelle A of genetic variant 1}
    \item{\code{allele.B1}}{Alelle B of genetic variant 1}
    \item{\code{allele.A2}}{Alelle A of genetic variant 2}
    \item{\code{allele.B2}}{Alelle B of genetic variant 2}
    \item{\code{env1.U}}{Exposure data for environment 1}
    \item{\code{env2.U}}{Exposure data for environment 2}
    \item{\code{int.U}}{Data for the interaction term }
  }
}

\examples{data(init.data)}

\keyword{datasets}
