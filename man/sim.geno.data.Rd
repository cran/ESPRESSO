\name{sim.geno.data}

\alias{sim.geno.data}

\title{Simulates genotypes for a genetic variant}

\description{
  Generates two alleles and combines them to form the genotype of a SNP under a binary or additive genetic model.
}

\usage{sim.geno.data(num.obs = 20000, MAF = 0.1, is.add = 0)}

\arguments{
  \item{num.obs}{Number of observations to simulate}
  \item{MAF}{Minor allele frequency of the variant}
  \item{is.add}{Genetic model of the variant}
}

\value{
  A dataframe containg the following variables:
  \item{allele.A}{Major allele}
  \item{allele.B}{Minor allele}
  \item{geno.U}{Genotype}
}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{sim.LDsnps}}, \code{\link{sim.LDgeno.data}}}

\examples{
# simulate genotypes for a binary SNP with a MAF of 0.1 
geno <- sim.geno.data(num.obs = 10000, 0.1, 0)
}

\keyword{ utilities }

