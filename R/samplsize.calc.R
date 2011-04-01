samplsize.calc <-
function(numcases=2000, numcontrols=8000, num.subjects=500, pheno.model=0, is.interaction=0, pval=1e-04, power=0.8, mean.model.z.geno, mean.model.z.env, mean.model.z.int)
{
   mean.model.z.geno1 <- mean.model.z.geno[1]
   mean.model.z.geno2 <- mean.model.z.geno[2]
   mean.model.z.env1 <- mean.model.z.env[1]
   mean.model.z.env2 <- mean.model.z.env[2]

   # CALCULATE Z STATISTIC THRESHOLD FOR DESIRED P-VALUE AND POWER
   z.pval <- qnorm(1-pval/2)
   z.power.required <- qnorm(power)+z.pval

   # ESTIMATE HOW MUCH THE SIMULATED STUDY SIZE NEEDS TO BE INFLATED OR SHRINKED 
   # IN ORDER TO OBTAIN A POWER OF 80%. THE RATIO OF Z STATISTIC REQUIRED FOR DESIRED 
   # POWER TO MEAN MODEL Z STATISTIC OBTAINED INDICATES RELATIVE CHANGE REQUIRED IN STANDARD ERROR. 
   # THIS CORRESPONDS TO RELATIVE CHANGE ON SCALE OF SQUARE ROOT OF SAMPLE SIZE. RATIO OF SAMPLE 
   # SIZE IS THEREFORE THIS RATIO SQUARED.
   sample.size.inflation.geno1 <- (z.power.required/mean.model.z.geno1)^2
   sample.size.inflation.geno2 <- (z.power.required/mean.model.z.geno2)^2

   # MULTIPLY THE ACTUAL NUMBER OF CASES AND CONTROLS BY THIS
   # SQUARED RATIO TO GET THE REQUIRED NUMBER OF CASES AND CONTROLS
   # FOR THE DESIRED POWER
   a <- round(numcases*sample.size.inflation.geno1,0)
   b <- round(numcontrols*sample.size.inflation.geno1,0)
   c <- round(numcases*sample.size.inflation.geno2,0)
   d <- round(numcontrols*sample.size.inflation.geno2,0)

   # NUMBER OF SUBJECTS REQUIRED FOR QTL - GENO
   q1 <- round(num.subjects*sample.size.inflation.geno1,0)
   q2 <- round(num.subjects*sample.size.inflation.geno2,0)

   # ENVIRONMENTAL DETERMINANT AS FOR GENOTYPE
   sample.size.inflation.env1 <- (z.power.required/mean.model.z.env1)^2
   e <- round(numcases*sample.size.inflation.env1,0)
   f <- round(numcontrols*sample.size.inflation.env1,0)
   sample.size.inflation.env2 <- (z.power.required/mean.model.z.env2)^2
   g <- round(numcases*sample.size.inflation.env2,0)
   h <- round(numcontrols*sample.size.inflation.env2,0)

   # NUMBER OF SUBJECTS REQUIRED FOR QTL - ENV
   q3 <- round(num.subjects*sample.size.inflation.env1,0)
   q4 <- round(num.subjects*sample.size.inflation.env2,0)

   if(is.interaction != 0)
   {
      # INTERACTION - IF PRESENT - AS FOR GENOTYPE AND ENVIRONMENT
      model.power.int <- pnorm(mean.model.z.int-z.pval)
      sample.size.inflation.int <- (z.power.required/mean.model.z.int)^2
      i <- round(numcases*sample.size.inflation.int,0)
      j <- round(numcontrols*sample.size.inflation.int,0)
      # NUMBER OF SUBJECTS REQUIRED FOR QTL
      q5 <- round(num.subjects*sample.size.inflation.int,0)
   }else{
      model.power.int <- NA
      sample.size.inflation.int <- NA
      i <- NA
      j <- NA
      q5 <- NA
   }

   if(pheno.model==0){
      # under binary outcome misclassification shrinks ORs toward the null, get estimated ORs
      return(list(numcases.required.geno1=a, numcontrols.required.geno1=b, numcases.required.geno2=c, numcontrols.required.geno2=d, 
      numcases.required.env1=e, numcontrols.required.env1=f, numcases.required.env2=g, numcontrols.required.env2=h, numcases.required.int=i,
      numcontrols.required.int=j))
   }else{
      return(list(numsubjects.required.geno1=q1, numsubjects.required.geno2=q2,numsubjects.required.env1=q3, numsubjects.required.env2=q4, 
      numsubjects.required.int=q5))
   }

}

