\name{make.obs.geno}

\alias{make.obs.geno}

\title{Adds some error to genotype data}

\description{
  Simulates errors and adds it to the true data to obtain observed data. The alleles simulated by \code{sim.geno.data} are randomly misclassified and used to form new genotypes that represent the observed genotypes. 
}

\usage{
make.obs.geno(allele1, allele2, error.1.0=0.05, error.0.1=0.05, 
is.add=0, MAF=0.1) 
}

\arguments{
  \item{allele1}{Allele A} 
  \item{allele2}{Allele B}
  \item{error.1.0}{1 to 0 misclassification rate}
  \item{error.0.1}{0 to 1 misclassification rate}
  \item{is.add}{Genetic model: Set to 0 for binary model and 1 for additive model}
  \item{MAF}{Minor Allele frequency}
}

\value{
  A dataframe containing the following columns
  \item{genotyp.U}{observed genotypes}
  \item{allele.A.orig}{true A alleles}
  \item{allele.A.new}{observed A alleles}
  \item{allele.B.orig}{true B alleles}
  \item{allele.B.new}{observed B alleles}
}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{sim.geno.data}}}

\examples{
# Example 1
#
# simulate genotypes for two independent genetic variants
geno <- sim.geno.data(10000,0.1,0)
allele.A <- geno$allele.A
allele.B <- geno$allele.B

# randomly misclassify the above simulated alleles and form observed genotypes
obs.geno <- make.obs.geno(allele.A, allele.B, 0.05, 0.05, 0, 0)
}

\keyword{utilities}

