\name{general.params}

\alias{general.params}

\docType{data}

\title{Main parameters of the simulations}

\description{
  A table of scenarios (rows) and parameters (columns).
}
\usage{data(general.params)}

\format{
  A data frame with 24 observations on the following 18 variables.
  \describe{
    \item{\code{scenario.id}}{The id of the scenario (each row stores parameters for one scenario}
    \item{\code{seed.val}}{Seed value}
    \item{\code{numsims}}{Number of runs for each simulation}
    \item{\code{numcases}}{Number of cases under binary outcome}
    \item{\code{numcontrols}}{Number of controls under binary outcome}
    \item{\code{num.subjects}}{Number of subjects under continuous outcome}
    \item{interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
    \item{\code{pheno.model}}{Type of the outcome; 0 for binary and 1 for continuous}
    \item{\code{disease.prev}}{Prevalence of the binary outcome}
    \item{\code{or.int}}{Odds ration of the interaction}
    \item{\code{int.efkt}}{Interaction effect}
    \item{\code{RR.5.95}}{Baseline odss ratio for subject on 95 percent population centile versus 5 percent centile. This parameter reflects the heterogeneity in disease risk arising from determinants that have not been measured or have not been included in the model. If this parameter is set to 10, the implication is that a high risk subject (someone at the upper 95 percent centile of population risk) is, all else being equal, at 10 times the odds of developing disease compared to someone else who is at low risk (individual at the lower 5 percent centile of population risk).}
    \item{\code{p.val}}{Cut-off p-value defining statistical significance}
    \item{\code{power}}{Desired power}
    \item{\code{sensitivity.pheno}}{Sensitivity of the assessment of binary outcome}
    \item{\code{specificity.pheno}}{Specificity of the assessment of binary outcome}
    \item{\code{reliability.pheno}}{Reliability of the assessment of continuous outcome}
    \item{\code{sim.sesp.geno.env}}{Tells if sensitivity and specificty values should be simulated for genotypic and environmental exposures assessment; set to 1 to simulate}
  }
}

\examples{data(general.params)}

\keyword{datasets}
