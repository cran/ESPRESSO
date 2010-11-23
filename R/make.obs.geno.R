make.obs.geno <-
function (allele1, allele2, error.1.0=0.05, error.0.1=0.05, is.add=0, MAF=0.1) 
{
   allele.A <- allele1
   allele.B <- allele2
   is.add.gene <- is.add
   mean.add <- (2*MAF*(1-MAF)+2*(MAF^2))
   mean.bin <- (2*MAF*(1-MAF)+(MAF^2))
   allele.A.orig <- allele.A
   allele.B.orig <- allele.B

   allele.A.new <- misclassify(allele.A, error.1.0, error.0.1)
   allele.B.new <- misclassify(allele.B, error.1.0, error.0.1)

   genotyp.new <- allele.A.new+allele.B.new

   if(is.add.gene==0){
      genotyp.U <- genotyp.new > 0
      genotyp.U <- genotyp.U - mean.bin
   }
   if(is.add.gene==1){
      genotyp.U <- genotyp.new
      genotyp.U <- genotyp.U - mean.add
   }

   # return a dataframe 
   data.frame <- data.frame(genotyp.U, allele.A.orig, allele.A.new, allele.B.orig, allele.B.new)
}

