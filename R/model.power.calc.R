model.power.calc <-
function(is.interaction=0, pval=1e-04, mean.model.z.geno, mean.model.z.env, mean.model.z.int)
{
   mean.model.z.geno1 = mean.model.z.geno[1]
   mean.model.z.geno2 = mean.model.z.geno[2]
   mean.model.z.env1 = mean.model.z.env[1]
   mean.model.z.env2 = mean.model.z.env[2]
   mean.model.z.int = mean.model.z.int

   # CALCULATE EMPIRICAL POWER OF MODEL FOR SIMULATED SAMPLE SIZE
   z.pval <- qnorm(1-pval/2)

   # GENOTYPE
   a <- pnorm(mean.model.z.geno1-z.pval)
   b <- pnorm(mean.model.z.geno2-z.pval)

   # ENVIRONMENT
   c <- pnorm(mean.model.z.env1-z.pval)
   d <- pnorm(mean.model.z.env2-z.pval)

   # INTERACTION
   if(is.interaction != 0)
   {
      e <- pnorm(mean.model.z.int-z.pval)
   }else{
   			e <- NA
   }

   return(list(model.power.geno1=a, model.power.geno2=b, model.power.env1=c, model.power.env2=d, model.power.int=e))
}

