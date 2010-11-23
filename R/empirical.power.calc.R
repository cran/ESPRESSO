empirical.power.calc <-
function(is.interaction=0, pval, z.geno1.results,z.geno2.results,z.env1.results, z.env2.results, z.int.results)
{
   # CALCULATE Z STATISTIC THRESHOLD FOR DESIRED P-VALUE 
   z.pval <- qnorm(1-pval/2)

   # THE PROPORTION OF SIMULATIONS IN WHICH THE Z STATISTIC FOR THE PARAMETER OF INTEREST 
   # EXCEEDS THE Z STATISTIC FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE

   a <- round(mean(z.geno1.results > z.pval),3)
   b <- round(mean(z.geno2.results > z.pval),3)
   c <- round(mean(z.env1.results > z.pval),3)
   d <- round(mean(z.env2.results > z.pval),3)

   if(is.interaction != 0)
   {
      e <- round(mean(z.int.results > z.pval),3)
   }else{
      e <- NA
   }
   return(list(empirical.power.geno1=a, empirical.power.geno2=b, empirical.power.env1=c, empirical.power.env2=d, empirical.power.int=e))
}

