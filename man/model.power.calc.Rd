\name{model.power.calc}

\alias{model.power.calc}

\title{Calculates theoretical power}

\description{Computes the power of the study from the set sample size and desired power.}

\usage{
model.power.calc(is.interaction = 0, pval = 1e-04, mean.model.z.geno, 
mean.model.z.env, mean.model.z.int)
}

\arguments{
  \item{is.interaction}{Type of interaction; 1 for gene-environment, 2 for gene-gene and 3 for environment-environment interaction}
  \item{pval}{Cut-off p-value defining statistical significance}
  \item{mean.model.z.geno}{mean z-statistics of the two genetic determinants}
  \item{mean.model.z.env}{mean z-statistics of the environmental determinants }
  \item{mean.model.z.int}{mean z-statistics of the interaction term}
}

\value{
  A list containing:
  \item{model.power.geno1}{Theoretical power under the main effect of genetic variant 1}
  \item{model.power.geno2}{Theoretical power under the main effect of genetic variant 2}
  \item{model.power.env1}{Theoretical power under the main effect environment 1}
  \item{model.power.env2}{Theoretical power under the main effect environment 2}
  \item{model.power.int}{Theoretical power under Interaction}
}

\author{Amadou Gaye, Paul Burton}

\seealso{\code{\link{empirical.power.calc}}}

\examples{

# sets the model
is.interaction <-  0 # no interaction

# cut-off p-value
pval <- 1e-04

# mean z-statistics for the effects of the two genetic variants
mean.model.z.geno <- c(0.269, 0.268)

# mean z-statistics for environmental exposures
mean.model.z.env <- c(2.508,  2.512)

# mean z-statistics for the interaction part
mean.model.z.int <- NA

# calculate the power reached under the initial sample size and desired power
power <- model.power.calc(is.interaction,pval,mean.model.z.geno, 
mean.model.z.env,mean.model.z.int)

}

\keyword{utilities}

