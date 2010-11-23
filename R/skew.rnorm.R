skew.rnorm <-
function (num.obs=20000, mean=0, sd=1, skewness=0) 
{
    distrib1 <- rnorm(num.obs)
    distrib2 <- rnorm(num.obs)
    compare <- (distrib2 > skewness * distrib1)
    distrib1[compare] <- (-distrib1[compare])
    y <- mean + sd * distrib1
    attr(y, "parameters") <- c(mean, sd, skewness)
    return(y)
}

