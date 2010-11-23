sim.interact.data <-
function (geno1.U, geno2.U, env1.U, env2.U, is.interaction=1) 
{ 					
   # CREATE GENE-ENVIRONMENT INTERACTION TERM 
   int1.U <- geno1.U*env1.U

   # CREATE GENE-GENE INTERACTION TERM 
   int2.U <- geno1.U*geno2.U

   # CREATE ENV-ENV INTERACTION TERM 
   int3.U <- env1.U*env2.U

   # CHOOSE PRODUCT INTERACTION TERM
   int.U <- 0
   if(is.interaction == 1){
      int.U <- int1.U
   }
   if(is.interaction == 2){
      int.U <- int2.U
   }
   if(is.interaction == 3){
      int.U <- int3.U
   }

   # return vector
   output <- int.U
}

