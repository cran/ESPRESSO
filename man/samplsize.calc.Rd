\name{samplsize.calc}

\alias{samplsize.calc}

\title{Calculates the sample size required to achieved the desired power}

\description{
   Estimates by how much the simulated study size needs to be inflated or shrinked in order to obtain the specified level of power. The ratio of z- statistic required for desired power to mean model z-statistic obtained indicates the relative change in standard error required. This corresponds to relative change on scale of square root of sample size. 
}

\usage{
samplsize.calc(numcases = 2000, numcontrols = 8000, num.subjects = 500, 
pheno.model = 0, is.interaction = 0, pval = 1e-04, power = 0.8, mean.model.z.geno, 
mean.model.z.env, mean.model.z.int)
}

\arguments{
  \item{numcases}{Number of cases when outcome is binary}
  \item{numcontrols}{Number of controls when outcome is binary}
  \item{num.subjects}{Number of subjects when outcome is continuous}
  \item{pheno.model}{Outcome type, 0 for binary and 1 for continuous}
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{pval}{Cut-off p-value defining statistical significance}
  \item{power}{Desired power}
  \item{mean.model.z.geno}{Ratio of mean \code{beta} estimate over mean \code{se} estimate for genetic main effect}
  \item{mean.model.z.env}{Ratio of mean \code{beta} estimate over mean \code{se} estimate for environment main effect}
  \item{mean.model.z.int}{Ratio of mean \code{beta} estimate over mean \code{se} estimate for Interaction effect}
}

\value{
  A list containing:
  \item{numcases.required.geno1}{Number of cases required to achieve the desired power under genetic variant 1 main effect}
  \item{numcontrols.required.geno1}{Number of controls required to achieve the desired power under genetic variant 1 main effect}
  \item{numcases.required.geno2}{Number of cases required to achieve the desired power under genetic variant 2 main effect}
  \item{numcontrols.required.geno2}{Number of controls required to achieve the desired power under genetic variant 2 main effect}
  \item{numcases.required.env1}{Number of cases required to achieve the desired power under environment 1 main effect}
  \item{numcontrols.required.env1}{Number of controls required to achieve the desired power under environment 1 main effect}
  \item{numcases.required.env2}{Number of cases required to achieve the desired power under environment 2 main effect}
  \item{numcontrols.required.env2}{Number of controls required to achieve the desired power under environment 2 main effect}
  \item{numcases.required.int}{Number of cases required to achieve the desired power under interaction effect}
  \item{numcontrols.required.int}{Number of controls required to achieve the desired power under interaction effect}
  \item{numsubjects.required.geno1}{Number of subjects required to achieve the desired power under genetic variant 1 main effect}
  \item{numsubjects.required.geno2}{Number of subjects required to achieve the desired power under genetic variant 2 main effect}
  \item{numsubjects.required.env1}{Number of subjects required to achieve the desired power under environment 1 main effect}
  \item{numsubjects.required.env2}{Number of subjects required to achieve the desired power under environment 2 main effect}
  \item{numsubjects.required.int}{Number of subjects required to achieve the desired power under interaction}
}

\seealso{\code{\link{empirical.power.calc}}, \code{\link{model.power.calc}}}

\examples{


# Example 1: sample size calculation for a binary outcome

# set outcome type to binary
pheno.model <-  0

# set the model
is.interaction <-  0 # no interaction

# cut-off p-value
pval <- 1e-04

# desired power
power <- 0.80

# mean z-statistics for genetic variants
mean.model.z.geno <- c(0.269, 0.268)

# mean z-statistics for environmental exposures
mean.model.z.env <- c(2.508, 2.512)

# mean z-statistics for the interaction part
mean.model.z.int <- NA

# calculate the sample size required under each model
sample.sizes.required <- samplsize.calc(2000,8000, NA, pheno.model, 
is.interaction,pval, power, mean.model.z.geno,mean.model.z.env,mean.model.z.int)

# Example 2: sample size calculation for a continuous outcome

# set outcome type to binary
pheno.model <-  1

# set the model
is.interaction <-  0 # no interaction

# cut-off p-value
pval <- 1e-04

# desired power
power <- 0.80

# mean z-statistics for genetic variants
mean.model.z.geno <- c(4.890, 4.872)

# mean z-statistics for environmental exposures
mean.model.z.env <- c(2.508, 2.512)

# mean z-statistics for the interaction part
mean.model.z.int <- NA

# calculate the sample size required under each model
sample.sizes.required <- samplsize.calc(NA,NA, 1000, pheno.model, 
is.interaction,pval, power, mean.model.z.geno,mean.model.z.env,mean.model.z.int)

}

\keyword{ utilities }
