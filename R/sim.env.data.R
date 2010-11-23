sim.env.data <-
function (num.obs=20000, env.expo=0, env.mean.lowlm=3.3, env.stdev.uplm=5, env.prev=0.1, skewness=0) 

{ 	
   # CREATE THE FIRST ENVIRONMENTAL COVARIATE 
   if(env.expo==0){   # BINARY DISTRIBUTION
      env.U <- rbinom(num.obs,1,env.prev)
      mean.env <- env.prev
      env.U <- env.U-mean.env
   }
   if(env.expo==1){   # NORMAL DISTRIBUTION (STDEV IS ALWAYS 1 - STANDARDIZATION)
      env.U <- skew.rnorm(num.obs, env.mean.lowlm, 1, skewness)
      env.U <- env.U-mean(env.U)     # mean centering
   }
   if(env.expo==2){  # UNIFORM DISTRIBUTION
			   if(env.mean.lowlm >= env.stdev.uplm)
			   {
			     cat("\n\nALERT!\n Uniform Distribution: The upper limit must be greater than the lower limit\n\n")
			   }
      env.U <- runif(num.obs, min=env.mean.lowlm, max=env.stdev.uplm)
      env.U <- env.U-mean(env.U)     # mean centering
   }

   # return a vector 
   env <- env.U
}

