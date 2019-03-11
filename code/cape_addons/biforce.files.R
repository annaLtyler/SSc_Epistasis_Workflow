#This function takes in a data.obj and makes biForce Toolbox files from it
#desired phenotypes and covariates should already be in place in the data.obj
#a pair of files is made for each phenotype, since missing values are not allowed
#in the phenotype file
#phenotypes are corrected for covariates
#genotypes must be in 0,1,2 format

biforce.files <- function(data.obj, scan.what = c("eigentraits", "raw.traits")){
	
	scan.check <- grep("e", scan.what)
	if(length(scan.check) > 0){
		 phenos <- data.obj$ET
		}else{
		phenos <- data.obj$pheno
		}
	
	#correct phenotypes for covariates
	num.pheno <- dim(phenos)[2]	
	geno <- data.obj$geno

	covar.locale <- which(data.obj$chromosome == 0)
	covar.mat <- data.obj$geno[,covar.locale,drop = FALSE]
	
	for(i in 1:num.pheno){
		pv <- phenos[,i]
		not.na.locale <- which(!is.na(pv))
		pv <- pv[not.na.locale]
		pv.geno <- geno[not.na.locale,]
		corrected.pv <- residuals(lm(phenos[,i]~covar.mat))
		sample.id <- paste("sample_", names(corrected.pv), sep = "")
		pheno.mat <- cbind(sample.id, corrected.pv)
		geno.mat <- geno[as.numeric(names(corrected.pv)),]
		geno.mat <- cbind(sample.id, geno.mat)
		colnames(geno.mat) <- c("sample_id", paste("snp_", 1:(dim(geno.mat)[2]-1), sep = ""))
		rownames(geno.mat) <- c(paste("sample_", 1:dim(geno.mat)[1], sep = ""))
		write.table(pheno.mat, file = paste("BiForce.", colnames(phenos)[i], ".Pheno.txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
		write.table(t(geno.mat), file = paste("BiForce.", colnames(phenos)[i], ".Geno.txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE)
		}
	
	
}