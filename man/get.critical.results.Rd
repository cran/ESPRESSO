\name{get.critical.results}

\alias{get.critical.results}

\title{Summarizes the main results}

\description{
  Gets the number of cases and controls or subjects and the empirical and theoretical power under each model and prints a summary on the screen.
}

\usage{
get.critical.results(scenario, is.interaction = 0, pheno.model = 0, 
is.add = c(0, 0), env.expo = c(0, 0), sample.sizes.required, empirical.power, 
model.power)
}

\arguments{
  \item{scenario}{Scenario number}
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{pheno.model}{Type of the outcome; 0 for binary and 1 for continuous}
  \item{is.add}{Genetic models of the two variants: 0 for binary model and 1 for additive model}
  \item{env.expo}{Models of the environmental exposures}
  \item{sample.sizes.required}{Number of cases and controls or number of subjects required to achieve the desired power}
  \item{empirical.power}{Estimated empirical power}
  \item{model.power}{Calculated theoretical power}
}

\value{
  A table containing the following variables:
  \item{models}{Model of each of the covariates}
  \item{numcases}{Number of cases required to achieve the desired power under each model}
  \item{numcases}{Number of controls required to achieve the desired power under each model}
  \item{powers1}{Estimated empirical power under each model}
  \item{powers1}{Power achieved under each model with the specified sample size}
  \item{models}{Model of each of the covariates}
  \item{numsubjects}{Number of subjects required to achieve the desired power under each model}
  \item{powers1}{Estimated empirical power under each model}
  \item{powers1}{Power achieved under each model with the specified sample size}

}

\author{Amadou Gaye}

\examples{

# scenario number
j <- 1

# main effect model
is.interaction <-  1

# outcome
pheno.model <- 1

# models of the genetic variants
is.add <- c(0,0)

# models of the environmental exposures
env.expo <- c(0,0)

# Estimated sample sizes required for a continuous outcome
sample.sizes.required <- c(1000, 1300, 2000, 2400, 7000)

# Estimated values for empirical power
empirical.power <- c(0.2,0.3,0.56,0.6,0.15)

# power values calculated from the set number of subjects
model.power <- c(0.18,0.27,0.58,0.59,0.17)

# return critical results and print a summary
res <- get.critical.results(j, is.interaction, pheno.model, is.add, env.expo, 
sample.sizes.required, empirical.power, model.power)
}

\keyword{utilities}

