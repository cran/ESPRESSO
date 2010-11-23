\name{sim.CC.data}

\alias{sim.CC.data}

\title{Simulates cases and controls}

\description{
  Generates affected and non affected subjects until the set sample size is achieved.
}

\usage{
sim.CC.data(num.obs = 20000, numcases = 2000, numcontrols = 8000, 
allowed.sample.size = 2e+07, is.interaction = 0, disease.prev = 0.1, 
MAF = c(0.1, 0.1), is.add = c(0, 0), R.target = 0.7, LD = 0, cov.mat.req, 
display = FALSE, or.geno = c(1.5, 1.5), env.expo = c(0, 0), 
env.mean.lowlm = c(3.3, 3.3), env.stdev.uplm = c(1, 1), env.prev = c(0.1, 0.1), 
or.env = c(1.5, 1.5), skewness = c(0, 0), or.int = 1.8, sigma.subject = 12.36, 
pheno.error = c(0, 0))
}

\arguments{
  \item{num.obs}{Number of observations to generate per iteration}
  \item{numcases}{Number of cases to simulate}
  \item{numcontrols}{Number of controls to simulate}
  \item{allowed.sample.size}{Maximum number of observations allowed}
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{disease.prev}{Prevalence of the binary outcome}
  \item{MAF}{Minor allele frequencies of the genetic variants}
  \item{is.add}{Genetic models of the variants: 0 for binary model and 1 for additive model}
  \item{R.target}{Correlation coefficient required if the alleles of the two genetic variants are in LD}
  \item{LD}{Sets independence or LD between the two genetic variants: 0 for independence and 1 for LD}
  \item{cov.mat.req}{The covariance matrix required to generate 2 genetic variants in LD}
  \item{display}{If TRUE, a summary is printed on screen}
  \item{or.geno}{Odds ratios of the genetic variants}
  \item{env.expo}{Models of the environmental exposures}
  \item{env.mean.lowlm}{Mean under quantitative-normal model and lower limit under quantitative-uniform model}
  \item{env.stdev.uplm}{Standard deviation under quantitative-normal model and upper limit under quantitative-uniform model}
  \item{env.prev}{Prevalences of the environmental determinants}
  \item{or.env}{Odds ratios of the environmental determinants}
  \item{skewness}{Determines skewness under quantitative-normal model; right-skewed distribution if set to a positive value and left-skewed when set to a negative value}
  \item{or.int}{Odds ration of the interaction}
  \item{sigma.subject}{Baseline odds ratio for subject on 95 percent population centile versus 5 percent centile. This parameter reflects the heterogeneity in disease risk arising from determinants that have not been measured or have not been included in the model}
  \item{pheno.error}{Phenotype misclassification rates}
}

\value{A matrix containing 11 variables}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{sim.geno.data}}, \code{\link{sim.LDgeno.data}}, \code{\link{sim.env.data}}, \code{\link{sim.subject.data}}, \code{\link{sim.pheno.bin}}}

\examples{

# number of cases 
numcases <- 2000

# number of controls
numcontrols <- 8000

# main effect model
is.interaction <-  0 

# generate cases and controls untill the set number of cases, 
# controls and sample size is achieved 
sim.matrix <- sim.CC.data(20000,numcases,numcontrols,20000000,is.interaction,0.1,
c(0.1,0.1),c(0,0),0.7,0,cov.mat.req=NULL, FALSE,c(1.5,1.5),c(0,0),c(3.3,3.3),c(1,1), 
c(0.1,0.1),c(1.5,1.5),c(0,0),1.8,12.36,c(0,0))

}

\keyword{utilities}

