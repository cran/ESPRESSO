sim.geno.data <-
function(num.obs=20000, MAF=0.1, is.add=0)
{
  is.add.gene <- is.add

  # CORRECTION TERM FOR MEAN CENTERING FOR ADDITIVE 
  mean.geno.add.gene<-(2*MAF*(1-MAF)+2*(MAF^2))

  # CORRECTION TERM FOR MEAN CENTERING FOR BINARY GENE
  mean.geno.binary<-(2*MAF*(1-MAF)+(MAF^2))

  # CREATE, CENTRE AND ROUND AN ADDITIVE GENETIC GENOTYPE COVARIATE WITH
  # APPROPRIATE MAF 
  allele.A <- rbinom(num.obs,1,MAF)
  allele.B <- rbinom(num.obs,1,MAF)

  # ACTUAL GENO IS SUM OF ALLELES IN ADDITIVE.GENETIC MODEL
  # AND IS 1 IF SUM OF ALLELES IS 1 OR GREATER IN THE BINARY MODEL 
  geno <- allele.A+allele.B

		if(is.add.gene==0)
		{
				geno.U <- geno > 0
				geno.U <- geno.U - mean.geno.binary
		}
		if(is.add.gene==1)
		{
				geno.U <- geno
				geno.U <- geno.U - mean.geno.add.gene
		}

  # return the data as a dataframe
  dtframe <- data.frame(allele.A, allele.B, geno.U)
}

