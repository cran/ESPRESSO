sim.geno.sesp <-
function(seed.val=333333, prevalence.exp=0.5, R2.target=0.8)
{
   if(seed.val>0){
     set.seed(seed.val)
   }

   # SET MAX NUMBER OF ITERATIONS AND CONVERGENCE CRITERION
   maxiter <- 7
   convergence.criterion <- 0.0005

   # START SMALL AND USE LARGER SAMPLE SIZES NEAR TO OPTIMIZATION
   sample.size.vector <- c(500000,500000,500000,1000000,2000000)

   # SET RELIABILITY VALUES (EQUAL TO NUMBER OF ITERATIONS)
   start.vector <- c(0.00001,0.75,0.9,0.95,0.99,0.999,0.99999)
   reliability <- start.vector

   # LENGTH OF THE RELIABILITY VECTOR AND MEDIAN OF VECTOR
   start.length <- length(start.vector)
   midpoint <- round((start.length/2)+1)

   # VECTORS TO RECORD R2 VALUES 
   R2 <- rep(NA,start.length)

   # MATRIX TO RECORD ESTIMATED RELIABILITIES
   store.reliability <- matrix(NA,nrow=start.length,ncol=(maxiter+1))

   # COLUMN 1 OF THE ABOVE MATRIX IS THE VECTOR OF SET RELIABILITIES VALUES
   store.reliability[,1] <- reliability

   # SET CONVERGENCE AND ITERATIONS VALUES TO START WITH
   convergence <- 1
   iter <- 0

   # LOOP UNTILL CONVERGENCE AND MAXIMUM NUMBER OF ITERATION ARE REACHED
   while((abs(convergence) > convergence.criterion) & (iter < maxiter)){

      # CHOOSE THE SAMPLE SIZE TO START WITH
      iter <- iter+1
      choose.sample.size <- iter
      if(iter > 5){choose.sample.size <- 5}
      sample.size <- sample.size.vector[choose.sample.size]

      # PRINT SOME INFO ON THE SCREEN
      cat("\n\nSTART OF INTERATION",iter,"\n")
      cat("Sample size this iteration =",sample.size,"\nSpeed slows with increasing sample size\n\n")

      # 
      for(j in 1:start.length){

         # PRINT INFO TO TELL WHERE WE ARE WITHIN ONE iTERATION
         print(c(iter,j))

         # SET RELIABILTY VALUE AND ERROR (1-RELIABITY)
         var.b <- reliability[j]
         var.w <- (1-var.b)
         
         # GENERATE TRUE AND OBESERVED PHENOTYPE DATA
         pheno.orig <- rnorm(sample.size,0,sqrt(var.b))
         error.obs1 <- rnorm(sample.size,0,sqrt(var.w))
         pheno.obs1 <- pheno.orig+error.obs1
         pheno.true <- pheno.orig

         # GET THRESHOLD FOR THE TRUE DATA (THE VALUE BELOW WHICH 'prevalence.exp'% OF THE POINTS FALL)
         affected.threshold.true <- quantile(pheno.true,prevalence.exp)

         # GET THRESHOLD FOR THE OBSERVED DATA 
         affected.threshold.obs1 <- quantile(pheno.obs1,prevalence.exp)

         # ASSIGNED 1 TO TRUE AFFFECTED AND 0 TO TRUE NON-AFFECTED
         affected.true <- (pheno.true < affected.threshold.true)*1

         # ASSIGNED 2 TO OBSERVED AFFFECTED AND 1 TO OBSERVED NON-AFFECTED
         affected.obs1 <- (pheno.obs1 < affected.threshold.obs1)+1

         # MAKE A TABLE THAT GIVES TRUE AND FALSE AFFECTED AND TRUE AND FALSE NON-AFFECTED
         tab.true.obs1 <- table(affected.true,affected.obs1)

         # GET THE TOTAL NUMBER OF AFFECTED (TRUE.POSITIVES + FALSE.NEGATIVES)
         total.affected <- tab.true.obs1[2,1] + tab.true.obs1[2,2]

         # GET THE TOTAL NUMBER OF 	NON-AFFECTED (TRUE.NEGATIVES + FALSE.POSITIVES)
         total.non.affected <- tab.true.obs1[1,1] + tab.true.obs1[1,2]

         # CALCULATE SENSITIVITY AND SPECIFICITY FOR THIS ITERATION
         sensitivity <- tab.true.obs1[2,2]/total.affected            
         specificity <- tab.true.obs1[1,1]/total.non.affected 

         # CALCULATE THE CORRELATION COEFICIENT BETWEEN TRUE AND OBSERVED AFFECTED
         R2[j] <- (cor.test(affected.true,affected.obs1)$estimate)^2

         # GET THE DATA FOR THE MEDIAN POINT 
         if(j==midpoint){
            pheno.orig.mid <- pheno.orig
            error.obs1.mid <- error.obs1
            pheno.obs1.mid <- pheno.obs1
            error.obs2.mid <- rnorm(sample.size,0,sqrt(var.w))
            pheno.obs2.mid <- pheno.orig+error.obs2.mid
            reliability.mid <- reliability[j]
            sensitivity.mid <- sensitivity
            specificity.mid <- specificity
            R2.mid <- R2[j]
         }
      }

      # ASSIGNED THE INITIAL SET RELIABILITY VALUES TO ANOTHER VECTOR
      reliability.save <- reliability

      # CHOOSE A GREATER RELIABILITY IF THE ESTIMATED R2 IS SMALLER THAN THE TARGET R2
      # AND CHOOSE A SMALLER RELIABILITY IF THE ESTIMATED R2 IS GREATER THAN THE TARGET R2
      for(t in 1:start.length){
         rev.t <- (start.length-t)+1
         # IF THE T.TH CALCULATED R2 < TARGET R2 REPLACE THE 1ST RELIABITY VALUE BY THE T.TH RELIABILITY.SAVE VALUE
         if(R2[t] < R2.target){ reliability[1] <- reliability.save[t] }
         # IF THE REVT.TH CALCULATED R2 > TARGET R2 REPLACE THE LAST RELIABITY VALUE BY THE REVT.TH RELIABILITY.SAVE VALUE
         if(R2[rev.t] > R2.target){ reliability[start.length] <- reliability.save[rev.t] }
      }
      cat("END OF ITERATION",iter,"\n")

      # CALCULATE THE AMOUNT BY WHICH THE 2ND T0 6TH RELIABILTY VALUES SHOULD BE INCREMENTED
      increment <- (reliability[start.length]-reliability[1])/(start.length-1)

      # CALCULATE CONVERGENCE VALUE AND PRINT IT ON SCREEN
      convergence <- abs(R2.mid-R2.target)
      cat("\n\n Convergence Value\n",convergence)

      # INCREMENT THE 2ND T0 6TH RELIABILTY VALUES
      for(k in 2:(start.length-1)){
         reliability[k] <- reliability[1]+(k-1)*increment
      }

      # RECORD THE RELIABITY IN THE MATRIX
      store.reliability[,(iter+1)] <- reliability
   }

   cov.table <- cov(cbind(pheno.obs1.mid,pheno.obs2.mid))
   ICC <- cov.table[1,2]
   cat("\n\nCovariance\n", ICC)

   # PRINT SUMMARY ON SCREEN
   cat("\n\n SUMMARY\n #######\n")
   cat("\n\n Random Seed\n",seed.val)
   cat("\n\n Sample size\n",sample.size)
   cat("\n\n Simulated reliability\n",reliability.mid)
   cat("\n\n Target R2\n",R2.target)
   cat("\n\n Empirical R2\n",R2.mid)
   cat("\n\n Empirical reliability\n",ICC)
   cat("\n\n Required      Required\n Sensitivity   Specificity\n",sensitivity.mid,"   ",specificity.mid,"\n\n")

   # RETURN THE SIMULATED SENSITIVITY AND SPECIFICITY
   return(c(sensitivity.mid, specificity.mid))
}

