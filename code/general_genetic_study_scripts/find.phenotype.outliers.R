#This script finds outliers in phenotypic distributions
#for further scrutiny. Given a table of phenotypes, this
#script finds the phenotypes that lie beyond 3 standard
#deviations and reports the individuals in each phenotype
#that are outliers.
#pheno.table is a matrix of phenotypes in which each column
#is a phenotype and each row is an individual
#Threshold is the number of standard deviations outside of
#which a value is considered an outlier. This defaults to 3
#label is a string that is appended to the name of the output
#table.

find.phenotype.outliers <- function(pheno.table, threshold = 3, label = NULL){
	
	find.outlier <- function(pheno){
		sd.pheno <- sd(as.numeric(pheno), na.rm = TRUE)
		mean.pheno <- mean(as.numeric(pheno), na.rm = TRUE)
		outlier.locale <- c(which(as.numeric(pheno) > (mean.pheno + sd.pheno)), which(as.numeric(pheno) < mean.pheno - sd.pheno))
		vals <- pheno[outlier.locale]
		return(vals)
		}
	
	all.outliers <- apply(pheno.table, 2, find.outlier)
	
}