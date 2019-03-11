#This script takes in the names of two 
#columns and adds a new phenotype that is
#either the difference or ratio between them
#it saves the new object in the given filename
#operation is the function you want to perform
#on the columns. It can be "-", "/", etc.

add.diff.or.ratio <- function(cross.obj, orig.pheno1, orig.pheno2, operation, new.name){
	
	#make sure the new phenotype isn't here already
	new.pheno.locale <- which(colnames(cross.obj$pheno) == new.name)
	if(length(new.pheno.locale) > 0){
		cat(new.name, "has already been added\n")
		return(cross.obj)
		}
	
	FUN <- match.fun(operation)
	results <- FUN(as.numeric(cross.obj$pheno[,orig.pheno1]), as.numeric(cross.obj$pheno[,orig.pheno2])) 
	
	new.pheno <- cbind(cross.obj$pheno, results)
	colnames(new.pheno) <- c(colnames(cross.obj$pheno), new.name)

	cross.obj$pheno <- new.pheno

	return(cross.obj)

	}
