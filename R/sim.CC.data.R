sim.CC.data <-
function(num.obs=20000, numcases=2000, numcontrols=8000, allowed.sample.size=20000000, is.interaction=0, disease.prev=0.1, 
MAF=c(0.1,0.1), is.add=c(0,0), R.target=0.7, LD=0, cov.mat.req, display=FALSE, or.geno=c(1.5,1.5), env.expo=c(0,0), env.mean.lowlm=c(3.3,3.3), env.stdev.uplm=c(1,1), env.prev=c(0.1,0.1), or.env=c(1.5,1.5), skewness=c(0,0), or.int=1.8, sigma.subject=12.36, pheno.error=c(0,0))
{
   # SET UP ZEROED COUNT VECTORS TO DETERMINE WHEN ENOUGH CASES AND CONTROLS HAVE BEEN GENERATED
   complete <- 0
   complete.absolute <- 0
   cases.complete <- 0
   controls.complete <- 0
   block <- 0

   # SET UP A MATRIX WITH ONE EMPTY ROW - DUMMY ROW
   sim.matrix <- matrix(ncol = 11)

   # SET LOOP COUNTER
   numloops <- 0

   while(complete==0 && complete.absolute==0)
     {
        
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
					   env1.U  <- sim.env.data(num.obs, env.expo[1], env.mean.lowlm[1], env.stdev.uplm[1], env.prev[1], skewness[1])
					   env2.U <- sim.env.data(num.obs, env.expo[2], env.mean.lowlm[2], env.stdev.uplm[2], env.prev[2], skewness[2])

					   # CREATE NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR
					   # WITH APPROPRIATE VARIANCE ON SCALE OF LOG-ODDS
					   subject.effect <- sim.subject.data(num.obs, sigma.subject)

					   # CREATE INTERACTION TERM 
					   int.U <- sim.interact.data(geno1.U, geno2.U, env1.U, env2.U, is.interaction) 

        # CREATE OUTCOME DATA
        pheno.data <- sim.pheno.bin(num.obs, is.interaction, disease.prev, geno1.U, geno2.U, env1.U, env2.U, int.U, subject.effect, or.geno, or.env, or.int, pheno.error)

        pheno.original <- pheno.data$pheno.original
        pheno.U <- pheno.data$pheno.U


        # CREATE OUTPUT MATRIX WITH ONE ROW FOR EACH SUBJECT SIMULATED
        sim.matrix.temp <- cbind(pheno.original,pheno.U,geno1.U,geno2.U,allele.A1,allele.B1,allele.A2,allele.B2,env1.U,env2.U,int.U)

        # FOR THE 1ST LOOP DELETE THE DUMMY ROW AND APPEND THE GENERATED DATA
        if(numloops < 1)
        {    
           sim.matrix <- rbind(sim.matrix[-1,], sim.matrix.temp)
        }else{
           sim.matrix <- rbind(sim.matrix, sim.matrix.temp)
        }

        # SELECT OUT CASES
        sim.matrix.cases <- sim.matrix[pheno.U==1,]

        # SELECT OUT CONTROLS
        sim.matrix.controls <- sim.matrix[pheno.U==0,]

        # COUNT CASES AND CONTROLS
        cases.simulated <- dim(sim.matrix.cases)[1]
        controls.simulated <- dim(sim.matrix.controls)[1]

        # TEST IF THERE ARE AT LEAST ENOUGH CASES ALREADY SIMULATED
        # IF THERE ARE, DEFINE THE CASE ELEMENT OF THE DATA MATRIX
        if(cases.simulated>=numcases)
        {
           sim.matrix.cases <- sim.matrix.cases[1:numcases,]
           cases.complete <- 1
        }

        # TEST IF THERE ARE AT LEAST ENOUGH CONTROLS ALREADY SIMULATED
        # IF THERE ARE, DEFINE THE CONTROL ELEMENT OF THE DATA MATRIX
        if(controls.simulated>=numcontrols)
        {
           sim.matrix.controls <- sim.matrix.controls[1:numcontrols,]
           controls.complete <- 1
        }

        # HAVE WE NOW GENERATED THE REQUIRED NUMBER OF CASES AND CONTROLS?
        complete <- cases.complete*controls.complete		

        # HAVE WE EXCEEDED THE TOTAL SAMPLE SIZE ALLOWED?
        complete.absolute <- (((block+1)*num.obs)>=allowed.sample.size)
        if(complete.absolute==1) {sample.size.excess <- 1}

        # INCREMENT LOOP COUNTER
        numloops <- numloops + 1
    }

   # STACK FINAL DATA MATRIX WITH CASES FIRST
   sim.matrix <- rbind(sim.matrix.cases,sim.matrix.controls)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # REMOVE THE ORIGINAL PHENOTYPE (TRUE PHENOTYPES - COLUMN 2) TO MAKE SURE ONLY THE OBSERVED PHENOTYPES ARE USED
   sim.matrix <- sim.matrix[,c(1,3:12)]
   colnames(sim.matrix) <- c("id", "cc.U", "geno1.U", "geno2.U", "allele.A1","allele.B1", "allele.A2","allele.B2", "env1.U", "env2.U", "int.U")

   # RETURN A DATAFRAME
   mm <- data.frame(sim.matrix)
}

