get.observed.data <-
function(is.interaction=0, true.data, geno.error=c(0.05,0.05), is.add=c(0,0), MAF=c(0.1,0.1), env.expo=c(0,0), env.prev=c(0.1,0.1), env.error=c(0.15,0.15), reliability.env = c(0.8,0.8))
{
    sim.df <- true.data
    geno.error.1.0 <- geno.error[1]
    geno.error.0.1 <- geno.error[2]
    is.add.gene1 <- is.add[1]
    is.add.gene2 <- is.add[2] 
    MAF.1 <- MAF[1]
    MAF.2 <- MAF[2]
    env1.expo <- env.expo[1]
    env2.expo <- env.expo[2]
    env1.prev <- env.prev[1]
    env2.prev <- env.prev[2]
    env.error.1.0 <- env.error[1]
    env.error.0.1 <- env.error[2]
    reliability.env1 <- reliability.env[1]
    reliability.env2 <- reliability.env[2]


	   # GENE1
	   genotyp1 <- sim.df$geno1.U
	   observed.geno1.data <- make.obs.geno(sim.df$allele.A1, sim.df$allele.B1, geno.error.1.0, geno.error.0.1, is.add.gene1, MAF.1)
	   sim.df$geno1.U <- observed.geno1.data$genotyp.U
	   genotyp.obs1 <- observed.geno1.data$genotyp.U
	   allele.A1.obs <- observed.geno1.data$allele.A.new
	   allele.B1.obs <- observed.geno1.data$allele.B.new
	   allele.A1.true <- observed.geno1.data$allele.A.orig
	   allele.B1.true <- observed.geno1.data$allele.B.orig


	   # GENE2
	   genotyp2 <- sim.df$geno2.U
	   observed.geno2.data <-make.obs.geno(sim.df$allele.A2, sim.df$allele.B2, geno.error.1.0, geno.error.0.1, is.add.gene2, MAF.2)
	   sim.df$geno2.U <- observed.geno2.data$genotyp.U
	   genotyp.obs2 <- observed.geno2.data$genotyp.U
	   allele.A2.obs <- observed.geno2.data$allele.A.new
	   allele.B2.obs <- observed.geno2.data$allele.B.new
	   allele.A2.true <- observed.geno2.data$allele.A.orig
	   allele.B2.true <- observed.geno2.data$allele.B.orig

	   # FIRST ENVIRONMENT 
	   env1.data <- make.obs.env(sim.df$env1.U, env1.expo, env1.prev, env.error, reliability.env1)
	   sim.df$env1.U <- env1.data$environ.new

	   # SECOND ENVIRONMENT
	   env2.data <- make.obs.env(sim.df$env2.U, env2.expo, env2.prev, env.error,  reliability.env2)
	   sim.df$env2.U <- env2.data$environ.new

	   # UPDATE THE INTERACTION DATA USING DETERMINANTS WITH ERROR (OBSERVED DATA)
	   if(is.interaction == 1){                     # gene-environment interaction
	      sim.df$int.U <- sim.df$geno1.U*sim.df$env1.U
	   }
	   if(is.interaction == 2){                     # gene-gene interaction
	      sim.df$int.U <- sim.df$geno1.U * sim.df$geno2.U
	   }
	   if(is.interaction == 3){                     # env-env interaction
	      sim.df$int.U <- sim.df$env1.U * sim.df$env2.U
	   }
    
    # RETURN
    return(list(sim.df=sim.df))
}

