\name{sim.LDsnps}

\alias{sim.LDsnps}

\title{Simulates alleles for two biallelic SNPs in Linkage Disequilibrium}

\description{
  Generates alleles data for pre-spcified alleles frequencies. The covariance matrix required to achieved the desired LD is computed and used to produce a random vector from a bivariate normal distribution.
}

\usage{
sim.LDsnps(num.obs, maf.snp1 = 0.1, maf.snp2 = 0.1, R.target = 0.7, cov.mat.req, 
display = FALSE)
}

\arguments{
  \item{num.obs}{Number of observations to simulate}
  \item{maf.snp1}{Minor allele frequency of the first snp}
  \item{maf.snp2}{Minor allele frequency of the second snp}
  \item{R.target}{Correlation coefficient, desired level of LD}
  \item{cov.mat.req}{Covariance matrix required required to achieved the desired LD}
  \item{display}{If TRUE, a summary is printed on screen}
}

\value{
  A dataframe of two variables where the rows represent haplotypes 
  \item{snp1.allele }{allele data for the first snp }
  \item{snp2.allele }{allele data for the second snp }
}

\references{Montana, G. 2005, \code{HapSim: a simulation tool for generating haplotype data with
pre-specified allele frequencies and LD coefficients.}, Bioinformatics, \bold{vol. 21 (23)}, pp.4309-4311.}

\author{Amadou Gaye}

\seealso{\code{\link{sim.LDgeno.data}}}

\examples{
# desired LD
R.target <- 0.8

# MAF of the first snp
maf.1 <- 0.1

# MAF of the second snp
maf.2 <- 0.1

# the covariance matrix required to achieve the desired LD
cor.mat <- matrix(c(1,R.target,R.target,1),2,2) # correlation mtrix
cov.mat.req <- make.cov.mat(cor.mat, c(1-maf.1, 1-maf.2))

# if the required covariance matrix is not positive-definite get 
# the nearest positive-definite matrix (tolerance = 1e-06)
if(!is.posdef(cov.mat.req, 0.000001)){
  cov.mat.req <- make.posdef(cov.mat.req, 0.000001)
}

# generate allele data 
alleles <- sim.LDsnps(10000, maf.1, maf.2, R.target, cov.mat.req, TRUE)
}

\keyword{utilities}

