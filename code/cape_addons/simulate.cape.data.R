#This function simulates normally distributed 
#phenotype data with a given correlation and
#N independent genotypes. It puts the data
#into a structure that can be processed by 
#cape pheno.cor is a single number specifying 
#the average correlation between the phenotypes.
#The function creates an intercross data set
#by default.

simulate.cape.data <- function(num.ind = 100, num.pheno = 2, pheno.cor = 0.7, num.geno = 100, genotypes = c(0, 0.5, 1), geno.balance = c(1,1,1)){
				
	pheno <- matrix(NA, nrow = num.ind, ncol = num.pheno)
	
	for(i in 1:num.pheno){
		if(i == 1){
			pheno[,i] <- rnorm(num.ind)
			# pheno[,i] <- rpois(num.ind, lambda = 4)
			pheno[,i] <- rexp(num.ind, rate = 1)
			}else{
			pheno[,i] <- pheno[,1] + (pheno.cor*rnorm(num.ind))
			# pheno[,i] <- pheno[,1] + (pheno.cor*rpois(num.ind, lambda = 4))
			pheno[,i] <- pheno[,1] + (pheno.cor*rexp(num.ind, rate = 1))
			}
		}
	colnames(pheno) <- paste("phenotype", 1:num.pheno, sep = "")
	
	new.geno <- NULL
	for(i in 1:length(genotypes)){
		new.geno <- c(new.geno, rep(genotypes[i], geno.balance[i]*100))
		}
	genotypes <- new.geno
	
	geno <- matrix(NA, nrow = num.ind, ncol = num.geno)
	for(i in 1:num.geno){
		geno[,i] <- sample(genotypes, num.ind, replace = TRUE)
		}	
	colnames(geno) <- marker.names <- 1:num.geno
	
	
	#combine the objects into a cape object
	
	data.obj <- list(pheno, geno, marker.names, rep(1, num.geno), 1:num.geno)
	names(data.obj) <- c("pheno", "geno", "marker.names", "chromosome", "marker.location")
	return(data.obj)
	
}