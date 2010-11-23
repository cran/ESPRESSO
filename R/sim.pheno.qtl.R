sim.pheno.qtl <-
function(num.subjects=10000, is.interaction=0, geno1.U, geno2.U, env1.U, env2.U, int.U, geno.efkt=c(0.25,0.25), env.efkt=c(0.25,0.25), int.efkt=0.5, reliability.pheno)
{  
   # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   alpha <- 0 

   beta.geno1 <-	geno.efkt[1]
			beta.geno2 <-	geno.efkt[2]
			beta.env1 <-	env.efkt[1]
			beta.env2 <-	env.efkt[2]
			beta.int <-	int.efkt
   num.obs <- num.subjects

   # GENERATE SIMULATED OUTCOMES UNDER THE CHOSEN MODEL
   # MAIN EFFECTS MODEL
   if(is.interaction == 0)
   {
      lp <- alpha+beta.geno1*geno1.U+beta.geno2*geno2.U+beta.env1*env1.U+beta.env2*env2.U
   }else{
   # MODEL WITH INTERACTION
      if(is.interaction==1)
      {
         lp <- alpha+beta.geno1*geno1.U+beta.env1*env1.U+beta.int*int.U  
      }
      if(is.interaction==2)
      {
         lp <- alpha+beta.geno1*geno1.U+beta.geno2*geno2.U+beta.int*int.U
      }
      if(is.interaction==3)
      {
         lp <- alpha+beta.env1*env1.U+beta.env2*env2.U+beta.int*int.U
      }
   }

   # GENERATE ORIGINAL PHENOTYPES (TRUE PHENOTYPE DATA) 
   pheno <- rnorm(num.obs,lp,1) # standardized SD (SD = 1)
   pheno.original <- pheno

   # USE THE RELIABITLITY OF PHENOTYPE ASSESSMENT TO COMPUTE THE VARIANCE OF MEASURED PHENOTYPES.
   # RELIABITY = (VAR.true.data/VAR.obs.data) AND VAR.obs.data = (VAR.true.data + VAR.measurement)
   # IT FOLLOWS THAT VAR.true.data + VAR.of.estimate = VAR.true.data/RELIABILITY AND THEN:
   # VAR.measurement = (VAR.true.data/RELIABILITY) - VAR.true.data
   var.m <- (1^2/reliability.pheno)-(1^2) # standardized SD (SD = 1)

   # ADD ERROR TO ORIGINAL PHENOTYPES TO GENERATE OBSERVED PHENOTYPE DATA
   pheno.U <- rnorm(num.obs,pheno.original,sqrt(var.m))

   # RETURN
   pheno.data <- data.frame(pheno.original, pheno.U)
}

