make.obs.env <-
function(env.data, env.expo=0, env.prev=0.1, env.error=c(0.15,0.15), reliability.env=0.8)
{ 
   environ.orig <- env.data
   misclass.rate.1.to.0 <- env.error[1]
   misclass.rate.0.to.1 <- env.error[2]
   numsubs <- length(environ.orig)

   if(env.expo == 0)
   { # BINARY EXPOSURE (THIS IS EQUIVALENT TO CASE-CONTROL MISCLASSIFICATION 
     # BUT OCCURS AFTER THE SUBJECTS HAVE BEEN SAMPLED INTO THE STUDY)

     # TURN ENV DATA TO "0s" AND "1s" BEFORE USING THE "MISCLASS" FUNCTION
     mean.env <- env.prev
     environ.temp <- env.data + mean.env
         
     # MISCLASSIFY
     environ <- misclassify(environ.temp, misclass.rate.1.to.0, misclass.rate.0.to.1)

     # TURN IT BACK TO ITS INITIAL FORMAT (mean centred)
     environ.new <- environ - mean.env

   }else{
     if(env.expo==1){ 
        var.error <- (1^2/reliability.env)-(1^2) # standardized SD (SD = 1)
        # ADD ERROR TO ORIGINAL ENVIRONMENTAL EXPOSURE TO GENERATE OBSERVED DATA
        environ.n <- rnorm(numsubs, environ.orig, sqrt(var.error))
     }
     if(env.expo==2){ # UNIFORM
        env.expo.error <- rnorm(numsubs, 0, sd(environ.orig))
        environ.n <- env.data + env.expo.error
     }
     environ.new <- environ.n
   }

   # RETURNED BY THE FUNCTION
   dt <- data.frame(environ.orig, environ.new)
}

