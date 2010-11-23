sim.QTL.data <-
function(num.subjects, is.interaction=0, MAF=c(0.1,0.1), is.add=c(0,0), R.target=0.7, LD=0, cov.mat.req, display=FALSE, geno.efkt=c(0.25,0.25), env.expo=c(0,0), env.mean.lowlm=c(3.3,3.3), env.stdev.uplm=c(1,1), env.prev=c(0.1,0.1), env.efkt=c(0.25,0.25), skewness=c(0,0), int.efkt=0.5, reliability.pheno=0.9)
{
    num.obs <- num.subjects
	   MAF.1 <- MAF[1]
	   MAF.2 <- MAF[2]
	   is.add.gene1 <- is.add[1]
	   is.add.gene2 <- is.add[2]
	   
    if(LD == 0)
    {
				   # GENERATE THE FIRST GENOTYPE DATA
				   geno1.data <- sim.geno.data(num.obs, MAF.1, is.add.gene1)
				   allele.A1 <- geno1.data$allele.A
				   allele.B1 <- geno1.data$allele.B
				   geno1.U <- geno1.data$geno.U

				   # GENERATE THE SECOND GENOTYPE DATA 
				   geno2.data <- sim.geno.data(num.obs, MAF.2, is.add.gene2)
				   allele.A2 <- geno2.data$allele.A
				   allele.B2 <- geno2.data$allele.B
				   geno2.U <- geno2.data$geno.U
    }else{
       LDgeno.data <- sim.LDgeno.data(num.obs, MAF, is.add, R.target, cov.mat.req, display)
       allele.A1 <- LDgeno.data$allele.A1
       allele.B1 <- LDgeno.data$allele.B1
       allele.A2 <- LDgeno.data$allele.A2
       allele.B2 <- LDgeno.data$allele.B2
       geno1.U <- LDgeno.data$geno1.U
       geno2.U <- LDgeno.data$geno2.U
    }

    # GENERATE ENVIRONMENTAL DETERMINANTS DATA
				env1.U <- sim.env.data(num.obs, env.expo[1], env.mean.lowlm[1], env.stdev.uplm[1], env.prev[1], skewness[1])
	   env2.U <- sim.env.data(num.obs, env.expo[2], env.mean.lowlm[2], env.stdev.uplm[2], env.prev[2], skewness[2])

	   # CREATE INTERACTION TERM 
	   int.U <- sim.interact.data(geno1.U, geno2.U, env1.U, env2.U, is.interaction) 

    # GENERATE OUTCOME DATA
    pheno.data <- sim.pheno.qtl(num.obs, is.interaction, geno1.U, geno2.U, env1.U, env2.U, int.U, geno.efkt, env.efkt, int.efkt, reliability.pheno)

    pheno.original <- pheno.data$pheno.original
    pheno.U <- pheno.data$pheno.U

    # CREATE OUTPUT MATRIX WITH ONE ROW FOR EACH SUBJECT SIMULATED
    sim.matrix <- cbind(pheno.original,pheno.U,geno1.U,geno2.U,allele.A1,allele.B1,allele.A2,allele.B2,env1.U,env2.U,int.U)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # REMOVE THE ORIGINAL PHENOTYPE (TRUE PHENOTYPES - COLUMN 2) TO MAKE SURE ONLY THE OBSERVED PHENOTYPES ARE USED
   sim.matrix <- sim.matrix[,c(1,3:12)]
   colnames(sim.matrix) <- c("id", "pheno.U", "geno1.U", "geno2.U", "allele.A1","allele.B1", "allele.A2","allele.B2", "env1.U", "env2.U", "int.U")

   # RETURN A DATAFRAME
   mm <- data.frame(sim.matrix)
}

