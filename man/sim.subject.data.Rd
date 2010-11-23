\name{sim.subject.data}

\alias{sim.subject.data}

\title{Simulates the individual effect related to heterogeneity in disease risk}

\description{
 The variation in baseline disease risk is assumed to be normally distributed on the logistic scale. If this parameter is set to 10, the implication is
that a 'high risk' subject (someone at the upper 95 percent centile of population risk) is, all else being equal, at 10 times the odds of developing disease compared to someone else who is at 'low risk' (at the lower 5 percent centile of population risk).
}

\usage{sim.subject.data(num.obs = 20000, sigma.subject = 12.36)}


\arguments{
  \item{num.obs}{Number of observations to simulate}
  \item{sigma.subject}{Baseline odds ratio for subject on 95 percent population centile versus 5 percent centile. This parameter reflects the heterogeneity in disease risk arising from determinants that have not been measured or have not been included in the model}
}

\value{
  A numerical vector
}

\author{Amadou Gaye, Paul Burton}

\examples{
# generate subject effect data with a baseline OR of 10
subject.effect <- sim.subject.data(20000, 10)

}

\keyword{utilities}

