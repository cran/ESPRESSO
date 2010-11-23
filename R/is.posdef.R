is.posdef <-
function(matrix, tolerance=0.000001)
{

  # IF ONE OR MORE EIGEN VALUES < TOLERANCE VALUE, MATRIX IS NOT POSITIVE DEFINITE
  negative.eigens <- length(which(eigen(matrix)$values < tolerance))
  if(negative.eigens > 0){FALSE}else{TRUE}
  
}

