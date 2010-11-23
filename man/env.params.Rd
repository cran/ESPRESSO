\name{env.params}

\alias{env.params}

\docType{data}

\title{Parameters to simulate environmental exposures data}

\description{
  A table of scenarios (rows) and parameters (columns).
}
\usage{data(env.params)}

\format{
  A data frame with 24 observations on the following 21 variables.
  \describe{
    \item{\code{scenario.id}}{Scenario number}
    \item{\code{model.env1}}{Models of the first environmental exposure: 0 for binary and 1 quantitative-normal and 2 for quantitative-uniform}
    \item{\code{model.env2}}{Models of the first environmental exposure: 0 for binary and 1 quantitative-normal and 2 for quantitative-uniform}
    \item{\code{env1.variance}}{Variance of the observed measures for environment 1}
    \item{\code{env2.variance}}{Variance of the observed measures for environment 2}
    \item{\code{env1.error.stdev}}{Standard deviation of the error in the observed measures of environment 1}
    \item{\code{env2.error.stdev}}{Standard deviation of the error in the observed measures of environment 2}
    \item{\code{env1.prev}}{Prevalence of environment 1}
    \item{\code{env2.prev}}{Prevalence of environment 2}
    \item{\code{env1.mean.lowlm}}{Mean measure for environment 1 under quantitative-normal model and lower limit under quantitative-uniform model}
    \item{\code{env2.mean.lowlm}}{Mean measure for environment 2 under quantitative-normal model and lower limit under quantitative-uniform model}
    \item{\code{env1.stdev.uplm}}{Standard deviation under quantitative-normal model and upper limit under quantitative-uniform model, environment 1}
    \item{\code{env2.stdev.uplm}}{Standard deviation under quantitative-normal model and upper limit under quantitative-uniform model, environment 2}
    \item{\code{or.env1}}{Odds ratio for environment 1}
    \item{\code{env1.efkt}}{Effect size for environment 1}
    \item{\code{or.env2}}{Odds ratio for environment 1}
    \item{\code{env2.efkt}}{Effect size for environment 2}
    \item{\code{sensitivity.env}}{Sensitivity of the assessment to environment 1}
    \item{\code{specificity.env}}{Specificity of the assessment to environment 2}
    \item{\code{skewness1}}{Determines skewness under quantitative-normal model for environment 1; right-skewed distribution if set to a positive value and left-skewed when set to a negative value}
    \item{\code{skewness2}}{Determines skewness under quantitative-normal model for environment 2; right-skewed distribution if set to a positive value and left-skewed when set to a negative value}
  }
}

\examples{data(env.params)}

\keyword{datasets}
