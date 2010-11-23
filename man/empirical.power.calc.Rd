\name{empirical.power.calc}

\alias{empirical.power.calc}

\title{Calculates empirical power}

\description{
    Determines the proportion of simulations in which the z-statistic for the parameter of interest exceeds the z-statistic for the desired level of statistical significance. 
}

\usage{
empirical.power.calc(is.interaction = 0, pval, z.geno1.results, z.geno2.results, 
z.env1.results, z.env2.results, z.int.results)
}

\arguments{
  \item{is.interaction}{Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3}
  \item{pval}{Cut-off p-value defining statistical significance}
  \item{z.geno1.results}{Vector of z-statistics for the main effect of genetic variant 1}
  \item{z.geno2.results}{Vector of z-statistics for the main effect of genetic variant 2}
  \item{z.env1.results}{Vector of z-statistics for the main effect of environment 1}
  \item{z.env2.results}{Vector of z-statistics for the main effect of environment 2}
  \item{z.int.results}{Vector of z-statistics for interaction effect}
}

\value{
  A list containing:
  \item{empirical.power.geno1}{Empirical power under the main effect of genetic variant 1}
  \item{empirical.power.geno2}{Empirical power under the main effect of genetic variant 2}
  \item{empirical.power.env1}{Empirical power under the main effect environment 1}
  \item{empirical.power.env2}{Empirical power under the main effect environment 2}
  \item{empirical.power.int}{Empirical power under Interaction}
}

\author{Amadou Gaye, Paul Burton}


\seealso{\code{\link{samplsize.calc}}, \code{\link{model.power.calc}}}

\examples{

# set outcome type to binary
pheno.model <-  0

# set the model
is.interaction <-  0 # main effect model

# number of runs 
numsims <- 100

# empty vectors to store the results of each run
# gene 
beta.geno1.results<-rep(NA,numsims)
se.geno1.results<-rep(NA,numsims)
z.geno1.results<-rep(NA,numsims)
beta.geno2.results<-rep(NA,numsims)
se.geno2.results<-rep(NA,numsims)
z.geno2.results<-rep(NA,numsims)
# environment
beta.env1.results<-rep(NA,numsims)
se.env1.results<-rep(NA,numsims)
z.env1.results<-rep(NA,numsims)
beta.env2.results<-rep(NA,numsims)
se.env2.results<-rep(NA,numsims)
z.env2.results<-rep(NA,numsims)
# interaction
beta.int.results<-rep(NA,numsims)
se.int.results<-rep(NA,numsims)
z.int.results<-rep(NA,numsims)

# tracer to detect exceeding max allowable sample size
sample.size.excess <- 0

# generate and analyse datasets one at a time 
for(s in 1:numsims) 
{
   if(pheno.model == 0){# under binary outcome
     # generate cases and controls untill the required number of cases, 
     # controls and sample size is achieved 
     sim.matrix <- sim.CC.data(20000,2000,8000,20000000,0,0.1,c(0.1,0.1),c(0,0),0.7, 
     0,cov.mat.req=NULL, FALSE,c(1.5,1.5),c(0,0),c(3.3,3.3),c(1,1), 
     c(0.1,0.1),c(1.5,1.5),c(0,0),1.8,12.36,c(0,0))

  }else{# under quantitative outcome model
     # generate the specified number of subjects
     sim.matrix <- sim.QTL.data(1000,0,c(0.1,0.1),c(0,0), 
     0.7,0,cov.mat.req=NULL,FALSE,c(0.25,0.25),c(0,0),c(3.3,3.3),c(1,1),c(0.1,0.1), 
     c(0.25,0.25),c(0,0),0.5,0.9)
  }

   # add appropriate errors to produce observed genotypes 
   observed.data <- get.observed.data(0, sim.matrix,c(0.95,0.95), 
   c(0,0),c(0.1,0.1),c(0,0),c(0.1,0.1),c(0.85,0.85),c(0.18,0.18),c(0.04,0.04))
   sim.df <- observed.data$sim.df

   # data analysis
   glm.estimates <- regr.analysis(is.interaction, pheno.model, sim.df)

   # genetic variants estimates 
   beta.geno1.results[s] <- glm.estimates[1]
   se.geno1.results[s] <- glm.estimates[2]
   z.geno1.results[s] <- glm.estimates[3]
   beta.geno2.results[s] <- glm.estimates[4]
   se.geno2.results[s] <- glm.estimates[5]
   z.geno2.results[s] <- glm.estimates[6]

   # environment estimates
   beta.env1.results[s] <- glm.estimates[7]
   se.env1.results[s] <- glm.estimates[8]
   z.env1.results[s] <- glm.estimates[9]
   beta.env2.results[s] <- glm.estimates[10]
   se.env2.results[s] <- glm.estimates[11]
   z.env2.results[s] <- glm.estimates[12]
   
   # interaction estimates
   beta.int.results[s] <- glm.estimates[13]
   se.int.results[s] <- glm.estimates[14]
		 z.int.results[s] <- glm.estimates[15]  

   # print tracer after every nth dataset created
   #if(s %% trace.interval ==0)cat("\n",s,"of",numsims,"completed in scenario",scenario.id)

}
cat("\n\n")

# summary of primary parameter estimates
# genetic variants 
mean.beta.geno1 <- round(mean(beta.geno1.results),3)
mean.se.geno1 <- round(sqrt(mean(se.geno1.results^2)),3)
mean.model.z.geno1 <- mean.beta.geno1/mean.se.geno1
mean.beta.geno2 <- round(mean(beta.geno2.results),3)
mean.se.geno2 <- round(sqrt(mean(se.geno2.results^2)),3)
mean.model.z.geno2 <- mean.beta.geno2/mean.se.geno2
mean.model.z.geno <- c(mean.beta.geno1/mean.se.geno1, mean.beta.geno2/mean.se.geno2)
# environments
mean.beta.env1 <- round(mean(beta.env1.results),3)
mean.se.env1 <-round(sqrt(mean(se.env1.results^2)),3)
mean.model.z.env1 <- mean.beta.env1/mean.se.env1
mean.beta.env2 <- round(mean(beta.env2.results),3)
mean.se.env2 <-round(sqrt(mean(se.env2.results^2)),3)
mean.model.z.env2 <- mean.beta.env2/mean.se.env2
mean.model.z.env <- c(mean.beta.env1/mean.se.env1, mean.beta.env2/mean.se.env2)
# interaction
if(is.interaction == 0){
   mean.beta.int <- NA
   mean.se.int <- NA
   mean.model.z.int <- NA
}else{
   mean.beta.int<-round(mean(beta.int.results),3)
   mean.se.int<-round(sqrt(mean(se.int.results^2)),3)
   mean.model.z.int <- mean.beta.int/mean.se.int
}

# calculate empirical power ie simply the proportion of simulations in which
# the z statistic for the parameter of interest exceeds the z statistic
# for the desired level of statistical significance
empirical.power <- empirical.power.calc(is.interaction,1e-04,z.geno1.results,
z.geno2.results,z.env1.results,z.env2.results, z.int.results)

}

\keyword{utilities}
