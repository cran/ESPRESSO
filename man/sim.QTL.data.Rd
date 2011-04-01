\name{sim.QTL.data}

\alias{sim.QTL.data}

\title{Simulates subjects for a continuous outcome}

\description{
  Generates the specified number of subjects using functions \code{\link{sim.geno.data}} or \code{\link{sim.LDgeno.data}}, \code{\link{sim.env.data}} and \code{\link{sim.pheno.bin}}
}


\usage{
sim.QTL.data(num.subjects, is.interaction = 0, MAF = c(0.1, 0.1), 
is.add = c(0, 0), R.target = 0.7, LD = 0, cov.mat.req, display = FALSE, 
geno.efkt = c(0.25, 0.25), env.expo = c(0, 0), env.mean.lowlm = c(3.3, 3.3), 
env.stdev.uplm = c(1, 1), env.prev = c(0.1, 0.1), env.efkt = c(0.25, 0.25), 
skewness = c(0, 0), int.efkt = 0.5, reliability.pheno=0.9)
}

\arguments{
  \item{num.subjects}{Number of subjects to simulate}
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{MAF}{Minor allele frequencies of the two genetic variants}
  \item{is.add}{Genetic models of the two variants: 0 for binary and 1 for additive}
  \item{R.target}{Correlation coefficient required if the alleles of the two genes are in Linkage Disequilibrium}
  \item{LD}{Is there an LD between the alleles of the two genes (0 for no LD and 1 for LD}
  \item{cov.mat.req}{Covariance matrix required to generate 2 genetic variants in LD}
  \item{display}{If TRUE, a summary is printed on screen}
  \item{geno.efkt}{Effects of the genetic variants}
  \item{env.expo}{Models of the environmental exposures}
  \item{env.mean.lowlm}{Mean under quantitative-normal model and lower limit under quantitative-uniform model}
  \item{env.stdev.uplm}{Standard deviation under quantitative-normal model and upper limit under quantitative-uniform model}
  \item{env.prev}{Prevalences of the environmental exposures}
  \item{env.efkt}{Effects of the environmental determinants}
  \item{skewness}{Determines skewness under quantitative-normal model; right-skewed distribution if set to a positive value and left-skewed when set to a negative value}
  \item{int.efkt}{Effect of the interaction}
  \item{reliability.pheno}{Reliability of the phenotype assessment}
}


\value{A matrix containg 11 variables}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{sim.geno.data}}, \code{\link{sim.LDgeno.data}}, \code{\link{sim.env.data}} and \code{\link{sim.pheno.qtl}}}

\examples{

# number of subjects
num.subjects <- 500

# main effect model
is.interaction <-  0 

# generate cases and controls untill the set number of cases, 
# controls and sample size is achieved 
sim.QTL.data(num.subjects, is.interaction, c(0.1, 0.1), c(0, 0), 0.7, 0, 
cov.mat.req=NULL, FALSE, c(0.25, 0.25), c(0, 0), c(3.3, 3.3), c(1, 1), c(0.1, 0.1), 
c(0.25, 0.25), c(0, 0), 0.5, 0.9)

}


\keyword{utilities}

