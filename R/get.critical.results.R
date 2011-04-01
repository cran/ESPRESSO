get.critical.results <-
function(scenario, is.interaction=0, pheno.model=0, is.add=c(0,0), env.expo=c(0,0), sample.sizes.required, empirical.power, model.power, mean.betas)
{
  is.add.geno1 <- is.add[1]
  is.add.geno2 <- is.add[2]
  env1.expo <- env.expo[1]
  env2.expo <- env.expo[2]
  ss.required <- sample.sizes.required
  emp.power <- empirical.power
  mod.power <- model.power

		if(is.add.geno1 == 0){mod1 <- "binary"}else{mod1 <- "additive"}
		if(is.add.geno2 == 0){mod2 <- "binary"}else{mod2 <- "additive"}
		if(env1.expo == 0){
		  mod3 <- "binary"
		}else
		{
		  if(env1.expo == 1){mod3 <- "normal.dist"}
		  if(env1.expo == 2){mod3 <- "uniform.dist"}
		}
		if(env2.expo == 0){
		  mod4 <- "binary"
		}else
		{
		  if(env2.expo == 1){mod4 <- "normal.dist"}
		  if(env2.expo == 2){mod4 <- "uniform.dist"}
		}
		if(is.interaction == 0){
		  mod5 <- "no"
		}else{
		  if(is.interaction == 1){mod5 <- "Gen1-Env1"}
		  if(is.interaction == 2){mod5 <- "Gen1-Gen2"}
		  if(is.interaction == 3){mod5 <- "Env1-Env2"}
		}
		models <- c(mod1,mod2,mod3,mod4,mod5)
  
  if(pheno.model == 0){
					numcases <- c()
					for(i in c(1,3,5,7,9)){
					  numcases <- append(numcases, ss.required[[i]])
					}

					numcontrols <- c()
					for(i in c(2,4,6,8,10)){
					  numcontrols <- append(numcontrols, ss.required[[i]])
					}
  }else{
					numsubjects <- c()
					for(i in 1:5){
					  numsubjects <- append(numsubjects, ss.required[[i]])
					}

  }

		powers1 <- c()
		for(i in 1:5){
		  x <- round(emp.power[[i]],2)
		  powers1 <- append(powers1, x)
		}

		powers2 <- c()
		for(i in 1:5){
		  y <- round(mod.power[[i]],2)
		  powers2 <- append(powers2, y)
		}

  if(pheno.model==0){
     # estimated ORs
     est.or.geno1 <- exp(mean.betas[1])
     est.or.geno2 <- exp(mean.betas[2])
     est.or.env1 <- exp(mean.betas[3])
     est.or.env2 <- exp(mean.betas[4])
     est.or.I <- exp(mean.betas[5])
     est.ORs <- c(est.or.geno1, est.or.geno2, est.or.env1, est.or.env2, est.or.I)

					cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
					cat("\nModels\n")
					cat("------\n")
					cat(" Outcome: binary; ")
					cat(" Gen1:",mod1,"; ")
					cat(" Gen2:",mod2,"; ")
					cat(" Env1:",mod3,"; ")
					cat(" Env2:",mod4,"; ")
					cat(" Interaction:",mod5)

					cat("\n\nNumber of cases required\n")
					cat("------------------------\n")
					cat(" Gen1:",numcases[1],"; ")
					cat(" Gen2:",numcases[2],"; ")
					cat(" Env1:",numcases[3],"; ")
					cat(" Env2:",numcases[4],"; ")
					cat(" Interaction:",numcases[5])

					cat("\n\nNumber of controls required\n")
					cat("---------------------------\n")
					cat(" Gen1:",numcontrols[1],"; ")
					cat(" Gen2:",numcontrols[2],"; ")
					cat(" Env1:",numcontrols[3],"; ")
					cat(" Env2:",numcontrols[4],"; ")
					cat(" Interaction:",numcontrols[5])

					cat("\n\nEmpirical power\n")
					cat("---------------\n")
					cat(" Gen1:",powers1[1],"; ")
					cat(" Gen2:",powers1[2],"; ")
					cat(" Env1:",powers1[3],"; ")
					cat(" Env2:",powers1[4],"; ")
					cat(" Interaction:",powers1[5])

					cat("\n\nModel power\n")
					cat("-----------\n")
					cat(" Gen1:",powers2[1],"; ")
					cat(" Gen2:",powers2[2],"; ")
					cat(" Env1:",powers2[3],"; ")
					cat(" Env2:",powers2[4],"; ")
					cat(" Interaction:",powers2[5])

					cat("\n\nEstimated ORs\n")
					cat("-----------\n")
					cat(" Gen1:",est.or.geno1,"; ")
					cat(" Gen2:",est.or.geno2,"; ")
					cat(" Env1:",est.or.env1,"; ")
					cat(" Env2:",est.or.env2,"; ")
					cat(" Interaction:",est.or.I)

					cat("\n\n---- END OF SUMMARY ----\n")

		  crit.res <- matrix(c(models,numcases,numcontrols,powers1,powers2,est.ORs),5,6)

  }else{
     # estimated ORs
     est.or.geno1 <- 'NA'
     est.or.geno2 <- 'NA'
     est.or.env1 <- 'NA'
     est.or.env2 <- 'NA'
     est.or.I <- 'NA'
     est.ORs <- c(est.or.geno1, est.or.geno2, est.or.env1, est.or.env2, est.or.I)

					cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
					cat("\nModels\n")
					cat("------\n")
					cat(" Outcome: quantitative; ")
					cat(" Gen1:",mod1,"; ")
					cat(" Gen2:",mod2,"; ")
					cat(" Env1:",mod3,"; ")
					cat(" Env2:",mod4,"; ")
					cat(" Interaction:",mod5)

					cat("\n\nNumber of subjects required\n")
					cat("------------------------\n")
					cat(" Gen1:",numsubjects[1],"; ")
					cat(" Gen2:",numsubjects[2],"; ")
					cat(" Env1:",numsubjects[3],"; ")
					cat(" Env2:",numsubjects[4],"; ")
					cat(" Interaction:",numsubjects[5])

					cat("\n\nEmpirical power\n")
					cat("---------------\n")
					cat(" Gen1:",powers1[1],"; ")
					cat(" Gen2:",powers1[2],"; ")
					cat(" Env1:",powers1[3],"; ")
					cat(" Env2:",powers1[4],"; ")
					cat(" Interaction:",powers1[5])

					cat("\n\nModel power\n")
					cat("-----------\n")
					cat(" Gen1:",powers2[1],"; ")
					cat(" Gen2:",powers2[2],"; ")
					cat(" Env1:",powers2[3],"; ")
					cat(" Env2:",powers2[4],"; ")
					cat(" Interaction:",powers2[5])

					cat("\n\nEstimated ORs\n")
					cat("-----------\n")
					cat(" Gen1:",est.or.geno1,"; ")
					cat(" Gen2:",est.or.geno2,"; ")
					cat(" Env1:",est.or.env1,"; ")
					cat(" Env2:",est.or.env2,"; ")
					cat(" Interaction:",est.or.I)

					cat("\n\n---- END OF SUMMARY ----\n")

		   crit.res <- matrix(c(models,numsubjects,powers1,powers2,est.ORs),5,5)
  }

		return(crit.res)
}

