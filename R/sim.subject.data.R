sim.subject.data <-
function (num.obs=20000, sigma.subject=12.36) 
{ 					
   # CONVERT BASELINE ODDS RATIO FROM 5th TO 95th PERCENTILES INTO THE
   # CORRESPONDING VARIANCE FOR A NORMALLY DISTRIBUTED RANDOM EFFECT 
   sigma2.subject.effect <- (log(sigma.subject)/(2*qnorm(0.95)))^2

   # CREATE NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR
   # WITH APPROPRIATE VARIANCE ON SCALE OF LOG-ODDS
   subject.effect <- rnorm(num.obs,0,sqrt(sigma2.subject.effect))

   # RETURN A VECTOR 
   output <- subject.effect
}

