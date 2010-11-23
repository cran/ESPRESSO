sim.LDgeno.data <-
function(num.obs=20000, MAF=c(0.1,0.1), is.add=c(0,0), R.target=0.7, cov.mat.req, display=FALSE)
{
     maf.snp1 <- MAF[1]
     maf.snp2 <- MAF[2]
     is.add.gene <- is.add

     # CORRECTION TERM FOR MEAN CENTERING FOR ADDITIVE AND BINARY GENE
     mean.geno.add.gene <- c(2*MAF[1]*(1-MAF[1])+2*(MAF[1]^2), 2*MAF[2]*(1-MAF[2])+2*(MAF[2]^2))
     mean.geno.binary <- c(2*MAF[1]*(1-MAF[1])+(MAF[1]^2), 2*MAF[2]*(1-MAF[2])+(MAF[2]^2))

     # GENERATE ALLELES
     A.alleles <- sim.LDsnps (num.obs, maf.snp1, maf.snp2, R.target, cov.mat.req, display)
     B.alleles <- sim.LDsnps (num.obs, maf.snp1, maf.snp2, R.target, cov.mat.req, display)

     # GENE 1
     allele.A1 <- A.alleles[,1]
     allele.B1 <- B.alleles[,1]
     geno1 <- allele.A1+allele.B1

     if(is.add.gene[1]==1)
     {
        geno1.U <- geno1
        geno1.U <- geno1.U-mean.geno.add.gene[1]
     }else{
        geno1.U <- geno1>0
        geno1.U <- geno1.U-mean.geno.binary[1]
     }

     # GENE 2
     allele.A2 <- A.alleles[,2]
     allele.B2 <- B.alleles[,2]     
     geno2 <- allele.A2+allele.B2
     if(is.add.gene[2]==1)
     {
        geno2.U <- geno2
        geno2.U <- geno2.U-mean.geno.add.gene[2]
     }else{
        geno2.U <- geno2 > 0
        geno2.U <- geno2.U-mean.geno.binary[2]
     }

   # RETURN A DATAFRAME
   dtframe <- data.frame(allele.A1, allele.B1, geno1.U, allele.A2, allele.B2, geno2.U)
 }

