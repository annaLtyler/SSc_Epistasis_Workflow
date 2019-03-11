#this script removes snps with low minor allele frequencies

remove.low.maf <- function(data.obj, maf = 0.05){
	
	geno <- data.obj$geno
	
	get.allele.freq <- function(genotypes){
		allele1 <- (length(which(genotypes == 0))*2 + length(which(genotypes == 0.5)))/(length(genotypes[!is.na(genotypes)])*2)
		allele2 <- (length(which(genotypes == 1))*2 + length(which(genotypes == 0.5)))/(length(genotypes[!is.na(genotypes)])*2)
		return(min(c(allele1, allele2)))
		}
	
	all.maf <- apply(geno, 2, get.allele.freq)
	low.maf <- which(all.maf <= maf)
	
	if(length(low.maf) > 0){
		data.obj$geno <- geno[,-low.maf]
		cat("The following markers were removed:\n")
		cat(colnames(geno)[low.maf], sep = "\n")
		}
	
	return(data.obj)
	
	
}