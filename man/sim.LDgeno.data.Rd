\name{sim.LDgeno.data}

\alias{sim.LDgeno.data}

\title{ Simulates genotypes for two genetic variants in LD}

\description{
  Generates alleles of two SNPs in LD and uses these alleles to form the genotypes of the genetic variants. Each variant can be binary or additive.
}

\usage{
sim.LDgeno.data(num.obs = 20000, MAF = c(0.1, 0.1), is.add = c(0, 0), R.target = 0.7, cov.mat.req, display = FALSE)
}

\arguments{
  \item{num.obs}{Number of observations to simulate}
  \item{MAF}{Minor allele frequencies of the two variants}
  \item{is.add}{Models of the two variants}
  \item{R.target}{Correlation coefficient, desired level of LD}
  \item{cov.mat.req}{Covariance matrix required required to achieved the desired LD}
  \item{display}{If TRUE, a summary is printed on screen}
}

\value{
  A dataframe containg the following variables:
  \item{allele.A1}{Major allele of variant 1}
  \item{allele.B1}{Minor allele of variant 1}
  \item{geno1.U}{Genotype of variant 1}
  \item{allele.A2}{Major allele of variant 2}
  \item{allele.B2}{Minor allele of variant 2}
  \item{geno2.U}{Genotype of variant 2}
}

\author{Amadou Gaye}

\seealso{\code{\link{sim.LDsnps}}, \code{\link{sim.geno.data}}}

\examples{

# desired LD
R.target <- 0.8

# MAFs of the two genetic variants
MAFs <- c(0.1,0.1)

# the covariance matrix required to generate 2 variants with the desired LD
cor.mat <- matrix(c(1,R.target,R.target,1),2,2) # correlation mtrix
cov.mat.req <- make.cov.mat(cor.mat, c(1-MAFs[1], 1-MAFs[2]))

# if the required covariance matrix is not positive-definite get 
# the nearest positive-definite matrix (tolerance = 1e-06)
if(!is.posdef(cov.mat.req, 0.000001)){
  cov.mat.req <- make.posdef(cov.mat.req, 0.000001)
}

# generate genotypes for two gentic variants in LD
LDgeno <- sim.LDgeno.data(10000, c(0.1, 0.1), c(0, 0), 0.8, cov.mat.req, TRUE)
}

\keyword{ utilities }

