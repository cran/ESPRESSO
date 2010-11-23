sim.pheno.bin <-
function(num.obs=20000, is.interaction=0, disease.prev=0.1, geno1.U, geno2.U, env1.U, env2.U, int.U, subject.effect, or.geno=c(1.5,1.5), or.env=c(1.5,1.5), or.int=1.8, pheno.error=c(0,0))
{  
   alpha <- log(disease.prev/(1-disease.prev))
   beta.geno1 <-	log(or.geno[1])
			beta.geno2 <-	log(or.geno[2])
			beta.env1 <-	log(or.env[1])
			beta.env2 <-	log(or.env[2])
			beta.int <-	log(or.int) 
   pheno.error.1.0 <- pheno.error[1]
   pheno.error.0.1 <- pheno.error[2]

   # GENERATE SIMULATED OUTCOMES UNDER THE CHOSEN MODEL
   # MAIN EFFECTS MODEL
   if(is.interaction == 0)
   {
       lp <- alpha+beta.geno1*geno1.U+beta.geno2*geno2.U+beta.env1*env1.U+beta.env2*env2.U+subject.effect
   }else{
   # MODEL WITH INTERACTION
      if(is.interaction==1)
      {
         lp <- alpha+beta.geno1*geno1.U+beta.env1*env1.U+beta.int*int.U+subject.effect
      }
      if(is.interaction==2)
      {
         lp <- alpha+beta.geno1*geno1.U+beta.geno2*geno2.U+beta.int*int.U+subject.effect
      }
      if(is.interaction==3)
      {
         lp <- alpha+beta.env1*env1.U+beta.env2*env2.U+beta.int*int.U+subject.effect
      }
   }

   # LOGISTIC TRANSFORMATION: MU IS THE PROBABILITY OF DISEASE
   mu <- exp(lp)/(1 + exp(lp))
   
   # GENERATE PHENOTYPE DATA
   pheno <- rbinom(num.obs,1,mu)

   # STORE THESE TRUE OUTCOMES (IE TRUE DISEASE-NON_DISEASE STATUS)
   pheno.original <- pheno

   # RANDOMLY MISSCLASSIFY TRUE DISEASE-NON_DISEASE STATUS USING SPECIFIED MISCLASSIFICATION RATES 
   pheno.U <- misclassify(pheno.original, pheno.error.0.1, pheno.error.1.0)

   # RETURN
   pheno.data <- data.frame(pheno.original, pheno.U)
}

