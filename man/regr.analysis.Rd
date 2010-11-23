\name{regr.analysis}

\alias{regr.analysis}

\title{carries out regression analysis}

\description{
  fits a conventional unconditional logistic regression model wit a binary or continuous phenotype as outcome and the genetic, enviromental, interaction determinants as covariates.
}

\usage{regr.analysis(is.interaction = 0, pheno.model = 0, sim.df)}

\arguments{
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{pheno.model}{Type of the outcome; 0 for binary and 1 for continuous}
  \item{sim.df}{A dataframe that contains covariates and outcome data}
}

\value{A vector containg the \code{beta}, \code{standard-error} and \code{z-statistic} of each of the covariates}

\author{Amadou Gaye, Paul Burton}

\examples{

# load a table containing covariates and binary outcome data
data(obs.data)

# binary outcome
pheno.model <-  0

# is there an interaction
is.interaction <- 0

# regression analysis
glm.estimates <- regr.analysis(is.interaction, pheno.model, obs.data)

}

\keyword{ utilities }
