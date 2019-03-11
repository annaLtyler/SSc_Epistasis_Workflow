#This script takes in a vector of phenotype 
#names and adds a new phenotype that is
#the log transform of those phenotypes

add.log.transform <- function(data.obj, pheno, new.name){
	

	for(i in 1:length(pheno)){
		
		#make sure the new phenotype isn't here already
		new.pheno.locale <- which(colnames(data.obj$pheno) == new.name[i])
		if(length(new.pheno.locale) > 0){
			cat(new.name, "has already been added\n")
			}else{
				pheno.locale <- get.col.num(data.obj$pheno, pheno[i])
				new.pheno <- log10(as.numeric(data.obj$pheno[,pheno.locale]))
				new.pheno.mat <- cbind(data.obj$pheno, new.pheno)
				colnames(new.pheno.mat) <- c(colnames(data.obj$pheno), new.name[i])
				data.obj$pheno <- new.pheno.mat
				}
		
			}
	
	return(data.obj)

	}
