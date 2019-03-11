#This function scales the genotypes matrix to be between 0 and 1

scale.geno <- function(data.obj, geno.obj){
	
	geno <- get.geno(data.obj, geno.obj)
	scaled.geno <- geno/max(geno, na.rm = TRUE)

	if(!missing(geno.obj)){
		geno.obj$geno <- scaled.geno
		return(geno.obj)
		}else{
		data.obj$geno <- scaled.geno	
		return(data.obj)
		}		

}