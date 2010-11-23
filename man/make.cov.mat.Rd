\name{make.cov.mat}

\alias{make.cov.mat}

\title{Generates the covariance matrix required to achieved the desired LD}

\description{
  Finds the covariance values required to achieve the specified frequency of major allele haplotype. 
}

\usage{make.cov.mat(cor.mat, freqs)}

\arguments{
  \item{cor.mat}{Correlation matrix}
  \item{freqs}{Major allele frequencies of the two snps}
}

\value{A 2X2 covariance matrix}

\references{Montana, G. 2005, \code{HapSim: a simulation tool for generating haplotype data with
pre-specified allele frequencies and LD coefficients.}, Bioinformatics, \bold{vol. 21 (23)}, pp.4309-4311.
}

\author{Amadou Gaye}

\seealso{\code{\link{sim.LDsnps}}}

\examples{

# MAF of the first snp
maf1 <- 0.1

# MAF of the second snp
maf2 <- 0.1

# target LD
R.target <- 0.8

# correlation matrix
cor.mat <- matrix(c(1,R.target,R.target,1),2,2) 

# covariance matrix required to generate 2 variants with the desired LD
cov.mat.req <- make.cov.mat(cor.mat, c(1-maf1, 1-maf2))

}

\keyword{utilities}

