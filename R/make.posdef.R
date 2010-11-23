make.posdef <-
function(matrix, tolerance=0.000001)
{
    # IF THE INPUT MATRIX IS NOT SQUARED STOP THE PROCESS 
    if(dim(matrix)[1] != dim(matrix)[2]){
      cat("\n Your matrix is not square! :\n")
      stop()
    }

    # COMPUTE THE EIGEN VALUES
    eig <- eigen(matrix)
    eigvals <- eig$values
 
    # COMPUTE THE NEAREST POSITIVE DEFINITE OF MATRIX, USING THE ALGORITHM OF NJ HIGHAM (1988)
    D = 2 * tolerance
    para.max = pmax(0, D - eigvals)
    pos.def.mat = eig$vectors %*% diag(para.max, dim(matrix)[1]) %*% t(eig$vectors)
    return(matrix + pos.def.mat)
}

