sim.env.sesp <-
function(seed.val=333333, prevalence.exp=0.5, reliability=0.8)
{
   # SET SEED
   set.seed(seed.val)

   # SET SAMPLE SIZE
   sample.size <- 4000000
   
   # GET RELIABILTY VALUE AND ERROR (1-RELIABITY)
   var.b <- reliability
   var.w <- (1-var.b)
         
   # GENERATE ONE VECTOR OF TRUE AND TWO VECTORS OBESERVED PHENOTYPE DATA
   pheno.orig <- rnorm(sample.size,0,sqrt(var.b))
   error.obs1 <- rnorm(sample.size,0,sqrt(var.w))
   error.obs2 <- rnorm(sample.size,0,sqrt(var.w))
   pheno.obs1 <- pheno.orig+error.obs1
   pheno.obs2 <- pheno.orig+error.obs2
   pheno.true <- pheno.orig

   # CALCULATE COVARIANCE BETWEEN THE 2 TWO VECTORS OF OBSERVED DATA
   cov.table <- cov(cbind(pheno.obs1,pheno.obs2))
   ICC <- cov.table[1,2]

   # GET THRESHOLD FOR THE TRUE DATA (THE VALUE BELOW WHICH 'prevalence.exp'% OF THE POINTS FALL)
   affected.threshold.true <- quantile(pheno.true,prevalence.exp)

   # GET THRESHOLD FOR THE OBSERVED DATA
   affected.threshold.obs1 <- quantile(pheno.obs1,prevalence.exp)
   
   # ASSIGNED 1 TO TRUE AFFFECTED AND 0 TO TRUE NON-AFFECTED
   affected.true <- (pheno.true<affected.threshold.true)*1

   # ASSIGNED 2 TO OBSERVED AFFFECTED AND 1 TO OBSERVED NON-AFFECTED
   affected.obs1 <- (pheno.obs1<affected.threshold.obs1)+1

   # MAKE A TABLE THAT GIVES TRUE AND FALSE AFFECTED AND TRUE AND FALSE NON-AFFECTED
   tab.true.obs1 <- table(affected.true,affected.obs1)

   # GET THE TOTAL NUMBER OF AFFECTED (TRUE.POSITIVES + FALSE.NEGATIVES)
   total.affected <- tab.true.obs1[2,1] + tab.true.obs1[2,2]

   # GET THE TOTAL NUMBER OF 	NON-AFFECTED (TRUE.NEGATIVES + FALSE.POSITIVES)
   total.non.affected <- tab.true.obs1[1,1] + tab.true.obs1[1,2]

   # CALCULATE SENSITIVITY AND SPECIFICITY 
   sensitivity <- tab.true.obs1[2,2]/total.affected            
   specificity <- tab.true.obs1[1,1]/total.non.affected 

   # PRINT SUMMARY ON SCREEN
   cat("\n\n SUMMARY\n #######\n")
   cat("\n Random Seed\n",seed.val)
   cat("\n\n Sample size\n",sample.size)
   cat("\n\n Simulated reliability\n",reliability)
   cat("\n\n Empirical reliability\n",ICC)
   cat("\n\n Required      Required\n Sensitivity   Specificity\n",sensitivity,"   ",specificity,"\n\n")

   # RETURN THE SIMULATED SENSITIVITY AND SPECIFICITY
   return(c(sensitivity, specificity))
}

