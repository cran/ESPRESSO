regr.analysis <-
function(is.interaction=0, pheno.model=0, sim.df)
{
	  # MAIN EFFECTS MODEL
	  if(is.interaction==0)
	  {
      # BINARY OUTCOME
      if(pheno.model == 0){
	        # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
	        mod.glm <- glm(cc.U~1+geno1.U+geno2.U+env1.U+env2.U,family=binomial,data=sim.df)
	        mod.sum <- summary(mod.glm)
      }
      # QUANTITATIVE OUTCOME     
      if(pheno.model == 1){
	        # FIT A GLM FOR QTL
	        mod.glm <- glm(pheno.U~1+geno1.U+geno2.U+env1.U+env2.U,family=gaussian, data=sim.df)
	        mod.sum <- summary(mod.glm)     
      }
	     beta.geno1 <- mod.sum$coefficients[2,1]
	     se.geno1 <- mod.sum$coefficients[2,2]
	     z.geno1 <- mod.sum$coefficients[2,3]

	     beta.geno2 <- mod.sum$coefficients[3,1]
	     se.geno2 <- mod.sum$coefficients[3,2]
	     z.geno2 <- mod.sum$coefficients[3,3]

	     beta.env1 <- mod.sum$coefficients[4,1]
	     se.env1 <- mod.sum$coefficients[4,2]
	     z.env1 <- mod.sum$coefficients[4,3]

	     beta.env2 <- mod.sum$coefficients[5,1]
	     se.env2 <- mod.sum$coefficients[5,2]
	     z.env2 <- mod.sum$coefficients[5,3]

	     # DUMMY TERMS FOR INTERACTION
	     beta.int <- NA
	     se.int <- NA
	     z.int <- NA
	  }else{
	  # MODEL WITH INTERACTION
      if(is.interaction==1) # GEN-ENV INTERACTION
      {
			      # BINARY OUTCOME
			      if(pheno.model == 0){
				        # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
				        mod.glm <- glm(cc.U~1+geno1.U+env1.U+int.U,family=binomial,data=sim.df)
				        mod.sum <- summary(mod.glm)
			      }
			      if(pheno.model == 1){
			      # QUANTITATIVE OUTCOME, FIT A GLM FOR QTL
				        mod.glm <- glm(pheno.U~1+geno1.U+env1.U+int.U,family=gaussian,data=sim.df)
				        mod.sum <- summary(mod.glm)     
         }    
         beta.geno1 <- mod.sum$coefficients[2,1]
         se.geno1 <- mod.sum$coefficients[2,2]
         z.geno1 <- mod.sum$coefficients[2,3]
         beta.env1 <- mod.sum$coefficients[3,1]
         se.env1 <- mod.sum$coefficients[3,2]
         z.env1 <- mod.sum$coefficients[3,3]
         beta.geno2 <- NA
         se.geno2 <- NA
         z.geno2 <- NA
         beta.env2 <- NA
         se.env2 <- NA
         z.env2 <- NA  
      }
      if(is.interaction==2) # GEN-GEN INTERACTION
      {
			      if(pheno.model == 0){
				        # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
				        mod.glm <- glm(cc.U~1+geno1.U+geno2.U+int.U,family=binomial,data=sim.df)
				        mod.sum <- summary(mod.glm)
			      }
			      # QUANTITATIVE OUTCOME
			      if(pheno.model == 1){
				        # FIT A GLM FOR QTL
				        mod.glm <- glm(pheno.U~1+geno1.U+geno2.U+int.U,family=gaussian,data=sim.df)
				        mod.sum <- summary(mod.glm)     
			      }
	        beta.geno1 <- mod.sum$coefficients[2,1]
	        se.geno1 <- mod.sum$coefficients[2,2]
	        z.geno1 <- mod.sum$coefficients[2,3]
	        beta.geno2 <- mod.sum$coefficients[3,1]
	        se.geno2 <- mod.sum$coefficients[3,2]
	        z.geno2 <- mod.sum$coefficients[3,3]   
	        beta.env1 <- NA
	        se.env1 <- NA
	        z.env1 <- NA
	        beta.env2 <- NA
	        se.env2 <- NA
	        z.env2 <- NA 
      }
      if(is.interaction==3) # ENV-ENV INTERACTION
      {
			      if(pheno.model == 0){
				        # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
				        mod.glm <- glm(cc.U~1+env1.U+env2.U+int.U,family=binomial,data=sim.df)
				        mod.sum <- summary(mod.glm)
			      }
			      # QUANTITATIVE OUTCOME
			      if(pheno.model == 1){
				        # FIT A GLM FOR QTL
				        mod.glm <- glm(pheno.U~1+env1.U+env2.U+int.U,family=gaussian,data=sim.df)
				        mod.sum <- summary(mod.glm)     
			      }
	        beta.env1 <- mod.sum$coefficients[2,1]
	        se.env1 <- mod.sum$coefficients[2,2]
	        z.env1 <- mod.sum$coefficients[2,3]
	        beta.env2 <- mod.sum$coefficients[3,1]
	        se.env2 <- mod.sum$coefficients[3,2]
	        z.env2 <- mod.sum$coefficients[3,3]
	        beta.geno1 <- NA
	        se.geno1 <- NA
	        z.geno1 <- NA
	        beta.geno2 <- NA
	        se.geno2 <- NA
	        z.geno2 <- NA 
      }
	     beta.int <- mod.sum$coefficients[4,1]
	     se.int <- mod.sum$coefficients[4,2]
	     z.int <- mod.sum$coefficients[4,3] 
 
	  }
   # RETURN A VECTOR
   return(c(beta.geno1,se.geno1,z.geno1,beta.geno2,se.geno2,z.geno2,beta.env1,se.env1,z.env1,beta.env2,se.env2,z.env2,beta.int,se.int,z.int))
}

