\name{obs.data}

\alias{obs.data}

\docType{data}

\title{Simulated genotypes, phenotypes and enviromental exposure data}

\description{
  A table of simulated \code{observed} data. The number of rows represents the number of cases and controls (if the outcome is binary) or the number of subjects (if the outcome is continuous).
}

\usage{data(obs.data)}
\format{
  A data frame with 10000 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{Scenario id}
    \item{\code{cc.U}}{Phenotypes}
    \item{\code{geno1.U}}{Genotypes for genetic variant 1}
    \item{\code{geno2.U}}{Genotypes for genetic variant 2}
    \item{\code{allele.A1}}{Alelle A of genetic variant 1}
    \item{\code{allele.B1}}{Alelle B of genetic variant 1}
    \item{\code{allele.A2}}{Alelle A of genetic variant 2}
    \item{\code{allele.B2}}{Alelle B of genetic variant 2}
    \item{\code{env1.U}}{Data for environment 1}
    \item{\code{env2.U}}{Data for environment 2}
    \item{\code{int.U}}{Data for the interaction term }
  }
}

\examples{data(obs.data)}

\keyword{datasets}
