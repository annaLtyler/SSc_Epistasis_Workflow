#This function replaces one phenotype name with another
#throughout the data object.

change.phenotype.name <- function(data.obj, orig.pheno, new.pheno){
	
	if(length(orig.pheno) != length(new.pheno)){
		stop("orig.pheno and new.pheno must be vectors of the same length.")
		}
	
	for(i in 1:length(data.obj)){
		for(j in 1:length(orig.pheno)){
			
			if(length(names(data.obj[[i]])) > 0){
				names.obj.locale <- which(names(data.obj[[i]]) == orig.pheno[j])
				if(length(names.obj.locale) > 0){
					names(data.obj[[i]])[names.obj.locale] <- new.pheno[j]
					}
				}
			
			for(d in 1:length(dimnames(data.obj[[i]]))){
				names.obj.locale <- which(dimnames(data.obj[[i]])[[d]] == orig.pheno[j])
				if(length(names.obj.locale) > 0){
					dimnames(data.obj[[i]])[[d]][names.obj.locale] <- new.pheno[j]
					}
				} #end looping through dimnames
			} #end looping through phenotype names
		
		
	} #end looping through data object elements
	
	return(data.obj)
	
	
}